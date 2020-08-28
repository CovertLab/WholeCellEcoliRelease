"""
TranscriptElongation

Transcription elongation sub-model.

TODO:
- use transcription units instead of single genes
- account for energy

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
from wholecell.utils import units

class TranscriptElongation(wholecell.processes.process.Process):
	""" TranscriptElongation """

	_name = "TranscriptElongation"

	# Constructor
	def __init__(self):
		super(TranscriptElongation, self).__init__()


	def initialize(self, sim, sim_data):
		super(TranscriptElongation, self).initialize(sim, sim_data)

		self.max_time_step = sim_data.process.transcription.max_time_step

		# Load parameters
		self.rnaPolymeraseElongationRateDict = sim_data.process.transcription.rnaPolymeraseElongationRateDict
		self.rnaIds = sim_data.process.transcription.rnaData['id']
		self.rnaLengths = sim_data.process.transcription.rnaData["length"].asNumber()
		self.rnaSequences = sim_data.process.transcription.transcriptionSequences
		self.ntWeights = sim_data.process.transcription.transcriptionMonomerWeights
		self.endWeight = sim_data.process.transcription.transcriptionEndWeight
		self.replichore_lengths = sim_data.process.replication.replichore_lengths
		self.chromosome_length = self.replichore_lengths.sum()

		# ID Groups of rRNAs
		self.idx_16S_rRNA = np.where(sim_data.process.transcription.rnaData['isRRna16S'])[0]
		self.idx_23S_rRNA = np.where(sim_data.process.transcription.rnaData['isRRna23S'])[0]
		self.idx_5S_rRNA = np.where(sim_data.process.transcription.rnaData['isRRna5S'])[0]

		# Mask for mRNAs
		self.is_mRNA = sim_data.process.transcription.rnaData['isMRna']

		# Views
		self.active_RNAPs = self.uniqueMoleculesView('active_RNAP')
		self.RNAs = self.uniqueMoleculesView('RNA')
		self.bulk_RNAs = self.bulkMoleculesView(self.rnaIds)
		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView(sim_data.moleculeIds.ppi)
		self.inactive_RNAPs = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.variable_elongation = sim._variable_elongation_transcription
		self.make_elongation_rates = sim_data.process.transcription.make_elongation_rates


	def calculateRequest(self):
		# Calculate elongation rate based on the current media
		current_media_id = self._external_states['Environment'].current_media_id

		self.rnapElongationRate = self.rnaPolymeraseElongationRateDict[current_media_id].asNumber(units.nt / units.s)

		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			self.rnapElongationRate,
			self.timeStepSec(),
			self.variable_elongation)

		# If there are no active RNA polymerases, return immediately
		if self.active_RNAPs.total_count() == 0:
			return

		# Determine total possible sequences of nucleotides that can be
		# transcribed in this time step for each partial transcript
		TU_indexes, transcript_lengths, is_full_transcript = self.RNAs.attrs(
			'TU_index', 'transcript_length', 'is_full_transcript')
		is_partial_transcript = np.logical_not(is_full_transcript)
		TU_indexes_partial = TU_indexes[is_partial_transcript]
		transcript_lengths_partial = transcript_lengths[is_partial_transcript]

		sequences = buildSequences(
			self.rnaSequences,
			TU_indexes_partial,
			transcript_lengths_partial,
			self.elongation_rates)

		sequenceComposition = np.bincount(
			sequences[sequences != polymerize.PAD_VALUE], minlength = 4)

		# Calculate if any nucleotides are limited and request up to the number
		# in the sequences or number available
		ntpsTotal = self.ntps.total_counts()
		maxFractionalReactionLimit = np.fmin(1, ntpsTotal / sequenceComposition)
		self.ntps.requestIs(maxFractionalReactionLimit * sequenceComposition)

		self.writeToListener(
			"GrowthLimits", "ntpPoolSize", self.ntps.total_counts())
		self.writeToListener(
			"GrowthLimits", "ntpRequestSize",
			maxFractionalReactionLimit * sequenceComposition)

		# Request full access to active RNAPs and RNAs
		self.active_RNAPs.request_access(self.EDIT_DELETE_ACCESS)
		self.RNAs.request_access(self.EDIT_DELETE_ACCESS)


	def evolveState(self):
		ntpCounts = self.ntps.counts()
		self.writeToListener("GrowthLimits", "ntpAllocated", ntpCounts)

		if self.active_RNAPs.total_count() == 0:
			return

		# Get attributes from existing RNAs
		TU_index_all_RNAs, length_all_RNAs, is_full_transcript, is_mRNA_all_RNAs, RNAP_index_all_RNAs = self.RNAs.attrs(
			'TU_index', 'transcript_length', 'is_full_transcript', 'is_mRNA',
			'RNAP_index')

		# Determine sequences of RNAs that should be elongated
		is_partial_transcript = np.logical_not(is_full_transcript)
		partial_transcript_indexes = np.where(is_partial_transcript)[0]
		TU_index_partial_RNAs = TU_index_all_RNAs[is_partial_transcript]
		length_partial_RNAs = length_all_RNAs[is_partial_transcript]
		is_mRNA_partial_RNAs = is_mRNA_all_RNAs[is_partial_transcript]
		RNAP_index_partial_RNAs = RNAP_index_all_RNAs[is_partial_transcript]

		sequences = buildSequences(
			self.rnaSequences,
			TU_index_partial_RNAs,
			length_partial_RNAs,
			self.elongation_rates)

		# Polymerize transcripts based on sequences and available nucleotides
		reactionLimit = ntpCounts.sum()
		result = polymerize(
			sequences,
			ntpCounts,
			reactionLimit,
			self.randomState,
			self.elongation_rates[TU_index_partial_RNAs])

		sequence_elongations = result.sequenceElongation
		ntps_used = result.monomerUsages

		# Calculate changes in mass associated with polymerization
		added_mass = computeMassIncrease(sequences, sequence_elongations,
			self.ntWeights)
		did_initialize = (length_partial_RNAs == 0) & (sequence_elongations > 0)
		added_mass[did_initialize] += self.endWeight

		# Calculate updated transcript lengths
		updated_transcript_lengths = length_partial_RNAs + sequence_elongations

		# Get attributes of active RNAPs
		coordinates, domain_index, direction, RNAP_unique_index = self.active_RNAPs.attrs(
			'coordinates', 'domain_index', 'direction', 'unique_index')

		# Active RNAP count should equal partial transcript count
		assert len(RNAP_unique_index) == len(RNAP_index_partial_RNAs)

		# All partial RNAs must be linked to an RNAP
		assert (np.count_nonzero(RNAP_index_partial_RNAs == -1) == 0)

		# Get mapping indexes between partial RNAs to RNAPs
		partial_RNA_to_RNAP_mapping, RNAP_to_partial_RNA_mapping = get_mapping_arrays(
			RNAP_index_partial_RNAs, RNAP_unique_index)

		# Rescale boolean array of directions to an array of 1's and -1's.
		# True is converted to 1, False is converted to -1.
		direction_rescaled = (2*(direction - 0.5)).astype(np.int64)

		# Compute the updated coordinates of RNAPs. Coordinates of RNAPs
		# moving in the positive direction are increased, whereas coordinates
		# of RNAPs moving in the negative direction are decreased.
		updated_coordinates = coordinates + np.multiply(
			direction_rescaled, sequence_elongations[partial_RNA_to_RNAP_mapping])

		# Reset coordinates of RNAPs that cross the boundaries between right
		# and left replichores
		updated_coordinates[
			updated_coordinates > self.replichore_lengths[0]
			] -= self.chromosome_length
		updated_coordinates[
			updated_coordinates < -self.replichore_lengths[1]
			] += self.chromosome_length

		# Update transcript lengths of RNAs and coordinates of RNAPs
		length_all_RNAs[is_partial_transcript] = updated_transcript_lengths
		self.RNAs.attrIs(transcript_length=length_all_RNAs)
		self.active_RNAPs.attrIs(coordinates=updated_coordinates)

		# Update added submasses of RNAs. Masses of partial mRNAs are counted
		# as mRNA mass as they are already functional, but the masses of other
		# types of partial RNAs are counted as nonspecific RNA mass.
		added_nsRNA_mass_all_RNAs = np.zeros_like(
			TU_index_all_RNAs, dtype=np.float64)
		added_mRNA_mass_all_RNAs = np.zeros_like(
			TU_index_all_RNAs, dtype=np.float64)

		added_nsRNA_mass_all_RNAs[is_partial_transcript] = np.multiply(
			added_mass, np.logical_not(is_mRNA_partial_RNAs))
		added_mRNA_mass_all_RNAs[is_partial_transcript] = np.multiply(
			added_mass, is_mRNA_partial_RNAs)

		self.RNAs.add_submass_by_name("nonspecific_RNA", added_nsRNA_mass_all_RNAs)
		self.RNAs.add_submass_by_name("mRNA", added_mRNA_mass_all_RNAs)

		# Determine if transcript has reached the end of the sequence
		terminal_lengths = self.rnaLengths[TU_index_partial_RNAs]
		did_terminate_mask = (updated_transcript_lengths == terminal_lengths)
		terminated_RNAs = np.bincount(
			TU_index_partial_RNAs[did_terminate_mask],
			minlength = self.rnaSequences.shape[0])

		# Assume transcription from all rRNA genes produce rRNAs from the first
		# operon. This is done to simplify the complexation reactions that
		# produce ribosomal subunits.
		n_total_16Srrna = terminated_RNAs[self.idx_16S_rRNA].sum()
		n_total_23Srrna = terminated_RNAs[self.idx_23S_rRNA].sum()
		n_total_5Srrna = terminated_RNAs[self.idx_5S_rRNA].sum()

		terminated_RNAs[self.idx_16S_rRNA] = 0
		terminated_RNAs[self.idx_23S_rRNA] = 0
		terminated_RNAs[self.idx_5S_rRNA] = 0

		terminated_RNAs[self.idx_16S_rRNA[0]] = n_total_16Srrna
		terminated_RNAs[self.idx_23S_rRNA[0]] = n_total_23Srrna
		terminated_RNAs[self.idx_5S_rRNA[0]] = n_total_5Srrna

		# Update is_full_transcript attribute of RNAs
		is_full_transcript_updated = is_full_transcript.copy()
		is_full_transcript_updated[
			partial_transcript_indexes[did_terminate_mask]] = True
		self.RNAs.attrIs(is_full_transcript=is_full_transcript_updated)

		# Remove partial transcripts that have finished transcription and are
		# not mRNAs from unique molecules (these are moved to bulk molecules)
		self.RNAs.delByIndexes(
			partial_transcript_indexes[np.logical_and(
				did_terminate_mask, np.logical_not(is_mRNA_partial_RNAs))])

		# Remove RNAPs that have finished transcription
		self.active_RNAPs.delByIndexes(
			np.where(did_terminate_mask[partial_RNA_to_RNAP_mapping]))

		n_terminated = did_terminate_mask.sum()
		n_initialized = did_initialize.sum()
		n_elongations = ntps_used.sum()

		# Get counts of new bulk RNAs
		n_new_bulk_RNAs = terminated_RNAs.copy()
		n_new_bulk_RNAs[self.is_mRNA] = 0

		# Update bulk molecule counts
		self.ntps.countsDec(ntps_used)
		self.bulk_RNAs.countsInc(n_new_bulk_RNAs)
		self.inactive_RNAPs.countInc(n_terminated)
		self.ppi.countInc(n_elongations - n_initialized)

		# Write outputs to listeners
		self.writeToListener(
			"TranscriptElongationListener", "countRnaSynthesized",
			terminated_RNAs)
		self.writeToListener(
			"TranscriptElongationListener", "countNTPsUSed", n_elongations)

		self.writeToListener("GrowthLimits", "ntpUsed", ntps_used)

		self.writeToListener(
			"RnapData", "actualElongations", sequence_elongations.sum())
		self.writeToListener("RnapData", "didTerminate", did_terminate_mask.sum())
		self.writeToListener(
			"RnapData", "terminationLoss",
			(terminal_lengths - length_partial_RNAs)[did_terminate_mask].sum())


	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		return inputTimeStep <= self.max_time_step


def get_mapping_arrays(x, y):
	"""
	Returns the array of indexes of each element of array x in array y, and
	vice versa. Assumes that the elements of x and y are unique, and
	set(x) == set(y).
	"""
	def argsort_unique(idx):
		"""
		Quicker argsort for arrays that are permutations of np.arange(n).
		"""
		n = idx.size
		argsort_idx = np.empty(n, dtype=np.int64)
		argsort_idx[idx] = np.arange(n)
		return argsort_idx

	x_argsort = np.argsort(x)
	y_argsort = np.argsort(y)

	x_to_y = x_argsort[argsort_unique(y_argsort)]
	y_to_x = y_argsort[argsort_unique(x_argsort)]

	return x_to_y, y_to_x
