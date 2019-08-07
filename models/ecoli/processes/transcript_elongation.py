#!/usr/bin/env python

"""
TranscriptElongation

Transcription elongation sub-model.

TODO:
- use transcription units instead of single genes
- account for energy

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

import numpy as np
from itertools import izip

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
from wholecell.utils import units
from wholecell.utils.random import stochasticRound
from wholecell.listeners.listener import WriteMethod

class TranscriptElongation(wholecell.processes.process.Process):
	""" TranscriptElongation """

	_name = "TranscriptElongation"

	# Constructor
	def __init__(self):
		super(TranscriptElongation, self).__init__()


	def initialize(self, sim, sim_data):
		super(TranscriptElongation, self).initialize(sim, sim_data)

		# Load parameters
		self.rnaPolymeraseElongationRateDict = sim_data.process.transcription.rnaPolymeraseElongationRateDict
		self.rnaIds = sim_data.process.transcription.rnaData['id']
		self.rnaLengths = sim_data.process.transcription.rnaData["length"].asNumber()
		self.rnaSequences = sim_data.process.transcription.transcriptionSequences
		self.ntWeights = sim_data.process.transcription.transcriptionMonomerWeights
		self.endWeight = sim_data.process.transcription.transcriptionEndWeight
		self.replichore_lengths = sim_data.process.replication.replichore_lengths
		self.chromosome_length = self.replichore_lengths.sum()

		# Get DNA polymerase elongation rate (used to determine mask for RNAPs
		# that are expected to collide with the replisome in the current
		# timestep)
		self.dnaPolyElngRate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s)))

		# ID Groups of rRNAs
		self.idx_16Srrna = np.where(sim_data.process.transcription.rnaData['isRRna16S'])[0]
		self.idx_23Srrna = np.where(sim_data.process.transcription.rnaData['isRRna23S'])[0]
		self.idx_5Srrna = np.where(sim_data.process.transcription.rnaData['isRRna5S'])[0]

		# Views
		self.active_RNAPs = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)
		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.active_replisomes = self.uniqueMoleculesView("active_replisome")
		self.fragmentBases = self.bulkMoleculesView(
			[id_ + "[c]" for id_ in sim_data.moleculeGroups.fragmentNT_IDs])


	def calculateRequest(self):
		# Calculate elongation rate based on the current media
		current_media_id = self._external_states['Environment'].current_media_id

		self.rnapElngRate = int(stochasticRound(self.randomState,
			self.rnaPolymeraseElongationRateDict[current_media_id].asNumber(units.nt / units.s) * self.timeStepSec()))

		# If there are no active RNA polymerases, return immediately
		if self.active_RNAPs.total_counts()[0] == 0:
			return

		# Determine total possible sequences of nucleotides that can be
		# transcribed in this time step for each polymerase
		TU_indexes, transcript_lengths = self.active_RNAPs.attrs(
			'TU_index', 'transcript_length')
		sequences = buildSequences(
			self.rnaSequences, TU_indexes, transcript_lengths,
			self.rnapElngRate)
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

		# Request full access to active RNAPs
		self.active_RNAPs.request_access(self.EDIT_DELETE_ACCESS)


	def evolveState(self):
		ntpCounts = self.ntps.counts()
		self.writeToListener("GrowthLimits", "ntpAllocated", ntpCounts)

		if self.active_RNAPs.total_counts()[0] == 0:
			return

		# Determine sequences that can be elongated
		TU_indexes, transcript_lengths, coordinates, domain_index, direction = self.active_RNAPs.attrs(
			'TU_index', 'transcript_length', 'coordinates', 'domain_index', 'direction')
		sequences = buildSequences(
			self.rnaSequences, TU_indexes, transcript_lengths,
			self.rnapElngRate)

		# Polymerize transcripts based on sequences and available nucleotides
		reactionLimit = ntpCounts.sum()
		result = polymerize(
			sequences, ntpCounts, reactionLimit, self.randomState)
		sequenceElongations = result.sequenceElongation
		ntpsUsed = result.monomerUsages

		# Calculate updated transcript lengths and coordinates of RNAPs
		updated_lengths = transcript_lengths + sequenceElongations

		# Convert boolean array of directions to an array of 1's and -1's.
		# True is converted to 1, False is converted to -1.
		direction_converted = (2*(direction - 0.5)).astype(np.int64)

		# Compute the updated coordinates of RNAPs. Coordinates of RNAPs
		# moving in the positive direction are increased, whereas coordinates
		# of RNAPs moving in the negative direction are decreased.
		updated_coordinates = coordinates + np.multiply(
			direction_converted, sequenceElongations)

		# If there are active replisomes, construct mask for RNAPs that are
		# expected to collide with replisomes in the current timestep. If the
		# sign of the differences between the updated coordinates of replisomes
		# and RNAPs are opposite to the sign of the differences between the
		# original coordinates, a collision occurs and the RNAP is knocked off
		# the chromosome.
		# TODO (ggsun): This assumes that replisomes elongate at maximum rates.
		# 	Ideally this should be done in the reconciler.
		all_collisions = np.zeros_like(coordinates, dtype=np.bool)
		headon_collisions = np.zeros_like(coordinates, dtype=np.bool)

		if self.active_replisomes.total_counts()[0] > 0:
			domain_index_replisome, right_replichore, coordinates_replisome = self.active_replisomes.attrs(
				"domain_index", "right_replichore", "coordinates")

			elongation_length = np.ceil(
				self.dnaPolyElngRate * self.timeStepSec())

			for rr, coord_rep, dmn_idx in izip(right_replichore,
					coordinates_replisome, domain_index_replisome):
				if rr:
					coordinates_mask = (
						np.multiply(coordinates - coord_rep,
							coordinates - (coord_rep + elongation_length)) < 0)
				else:
					coordinates_mask = (
						np.multiply(coordinates - coord_rep,
							coordinates - (coord_rep - elongation_length)) < 0)

				all_collisions_mask = np.logical_and(
					domain_index == dmn_idx, coordinates_mask)
				all_collisions[all_collisions_mask] = True

				# Collisions are head-on if replisomes and RNAPs are going in
				# opposite directions, hence the exclusive OR
				headon_collisions_mask = np.logical_and(
					all_collisions_mask, np.logical_xor(direction, rr))
				headon_collisions[headon_collisions_mask] = True

		# Remaining collisions are codirectional
		codirectional_collisions = np.logical_and(
			all_collisions, ~headon_collisions)

		n_total_collisions = all_collisions.sum()
		n_headon_collisions = headon_collisions.sum()
		n_codirectional_collisions = codirectional_collisions.sum()

		# Get coordinates for where the collisions occur
		headon_collision_coordinates = coordinates[headon_collisions]
		codirectional_collision_coordinates = coordinates[
			codirectional_collisions]

		# Reset coordinates of RNAPs that cross the boundaries between right
		# and left replichores
		updated_coordinates[
			updated_coordinates > self.replichore_lengths[0]
			] -= self.chromosome_length
		updated_coordinates[
			updated_coordinates < -self.replichore_lengths[1]
			] += self.chromosome_length

		# Calculate changes in mass associated with polymerization
		added_rna_mass = computeMassIncrease(
			sequences, sequenceElongations, self.ntWeights)
		didInitialize = (transcript_lengths == 0) & (sequenceElongations > 0)
		added_rna_mass[didInitialize] += self.endWeight

		# Update attributes and submasses of active RNAPs
		self.active_RNAPs.attrIs(
			transcript_length=updated_lengths,
			coordinates=updated_coordinates)
		self.active_RNAPs.add_submass_by_name("RNA", added_rna_mass)

		if n_total_collisions > 0:
			# Remove polymerases that are projected to collide with replisomes
			self.active_RNAPs.delByIndexes(np.where(all_collisions)[0])

			# Increment counts of inactive RNA polymerases
			self.inactiveRnaPolys.countInc(n_total_collisions)

			# Get lengths of transcripts that were terminated prematurely as a
			# result of the collision
			incomplete_sequence_lengths = updated_lengths[all_collisions]

			# Increment counts of bases in incomplete transcripts
			incomplete_sequences = buildSequences(
				self.rnaSequences, TU_indexes[all_collisions],
				np.zeros(n_total_collisions, dtype=np.int64),
				incomplete_sequence_lengths.max())

			base_counts = np.zeros(4, dtype=np.int64)

			for sl, seq in izip(incomplete_sequence_lengths, incomplete_sequences):
				base_counts += np.bincount(seq[:sl], minlength = 4)

			self.fragmentBases.countsInc(base_counts)
			self.ppi.countInc(n_total_collisions)

		# Determine if transcript has reached the end of the sequence
		terminalLengths = self.rnaLengths[TU_indexes]
		didTerminate = np.logical_and(
			updated_lengths == terminalLengths, ~all_collisions)
		terminatedRnas = np.bincount(
			TU_indexes[didTerminate], minlength = self.rnaSequences.shape[0])

		# Assume transcription from all rRNA genes produce rRNAs from the first
		# operon. This is done to simplify the complexation reactions that
		# produce ribosomal subunits.
		n_total_16Srrna = terminatedRnas[self.idx_16Srrna].sum()
		n_total_23Srrna = terminatedRnas[self.idx_23Srrna].sum()
		n_total_5Srrna = terminatedRnas[self.idx_5Srrna].sum()

		terminatedRnas[self.idx_16Srrna] = 0
		terminatedRnas[self.idx_23Srrna] = 0
		terminatedRnas[self.idx_5Srrna] = 0

		terminatedRnas[self.idx_16Srrna[0]] = n_total_16Srrna
		terminatedRnas[self.idx_23Srrna[0]] = n_total_23Srrna
		terminatedRnas[self.idx_5Srrna[0]] = n_total_5Srrna

		# Remove polymerases that have finished transcription from unique
		# molecules
		self.active_RNAPs.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = ntpsUsed.sum()

		# Update bulk molecule counts
		self.ntps.countsDec(ntpsUsed)
		self.bulkRnas.countsInc(terminatedRnas)
		self.inactiveRnaPolys.countInc(nTerminated)
		self.ppi.countInc(nElongations - nInitialized)

		# Write outputs to listeners
		self.writeToListener(
			"TranscriptElongationListener", "countRnaSynthesized",
			terminatedRnas)
		self.writeToListener(
			"TranscriptElongationListener", "countNTPsUSed", nElongations)

		self.writeToListener("GrowthLimits", "ntpUsed", ntpsUsed)

		self.writeToListener(
			"RnapData", "actualElongations", sequenceElongations.sum())
		self.writeToListener("RnapData", "didTerminate", didTerminate.sum())
		self.writeToListener(
			"RnapData", "terminationLoss",
			(terminalLengths - transcript_lengths)[didTerminate].sum())

		self.writeToListener(
			"RnapData", "n_total_collisions", n_total_collisions)
		self.writeToListener(
			"RnapData", "n_headon_collisions", n_headon_collisions)
		self.writeToListener(
			"RnapData", "n_codirectional_collisions",
			n_codirectional_collisions)
		self.writeToListener(
			"RnapData", "headon_collision_coordinates",
			headon_collision_coordinates, writeMethod=WriteMethod.fill)
		self.writeToListener(
			"RnapData", "codirectional_collision_coordinates",
			codirectional_collision_coordinates, writeMethod=WriteMethod.fill)
