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

		# Views
		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)
		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.active_replisomes = self.uniqueMoleculesView("active_replisome")
		self.fragmentBases = self.bulkMoleculesView(
			[id_ + "[c]" for id_ in sim_data.moleculeGroups.fragmentNT_IDs])


	def calculateRequest(self):
		# Calculate elongation rate based on the current nutrients
		current_nutrients = self._external_states['Environment'].nutrients

		self.rnapElngRate = int(stochasticRound(self.randomState,
			self.rnaPolymeraseElongationRateDict[current_nutrients].asNumber(units.nt / units.s) * self.timeStepSec()))

		# If there are no active RNA polymerases, return immediately
		if self.activeRnaPolys.total_counts()[0] == 0:
			return

		# Determine total possible sequences of nucleotides that can be
		# transcribed in this time step for each polymerase
		activeRnaPolys = self.activeRnaPolys.molecules_read_only()
		TU_indexes, transcript_lengths = activeRnaPolys.attrs(
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


	def evolveState(self):
		ntpCounts = self.ntps.counts()
		self.writeToListener("GrowthLimits", "ntpAllocated", ntpCounts)

		activeRnaPolys = self.activeRnaPolys.molecules()
		if len(activeRnaPolys) == 0:
			return

		# Determine sequences that can be elongated
		TU_indexes, transcript_lengths, coordinates, domain_index, direction = activeRnaPolys.attrs(
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
		activeRnaPolys.attrIs(
			transcript_length=updated_lengths,
			coordinates=updated_coordinates)
		activeRnaPolys.add_submass_by_name("RNA", added_rna_mass)

		# Get attributes of replisomes
		replisomes = self.active_replisomes.molecules_read_only()

		# If there are active replisomes, construct mask for RNAPs that are
		# expected to collide with replisomes in the current timestep. We
		# assume that if either the starting coordinate or the final coordinate
		# of the RNAP lies within the expected trajectory of replisomes, the
		# RNAP will collide with the replisome and fall off the chromosome.
		# TODO (ggsun): This assumes that replisomes elongate at maximum rates.
		# 	Ideally this should be done in the reconciler.
		collision_mask = np.zeros_like(coordinates, dtype=np.bool)

		if len(replisomes) > 0:
			domain_index_replisome, right_replichore, coordinates_replisome = replisomes.attrs(
				"domain_index", "right_replichore", "coordinates")

			elongation_length = np.ceil(
				self.dnaPolyElngRate * self.timeStepSec())

			for rr, coord, dmn_idx in izip(right_replichore,
					coordinates_replisome, domain_index_replisome):
				if rr:
					start_mask = np.logical_and(
						coordinates >= coord,
						coordinates <= coord + elongation_length)
					final_mask = np.logical_and(
						updated_coordinates >= coord,
						updated_coordinates <= coord + elongation_length)
				else:
					start_mask = np.logical_and(
						coordinates <= coord,
						coordinates >= coord - elongation_length)
					final_mask = np.logical_and(
						updated_coordinates <= coord,
						updated_coordinates >= coord - elongation_length)

				mask = np.logical_and(domain_index == dmn_idx,
					np.logical_or(start_mask, final_mask))
				collision_mask[mask] = True

		n_collisions = collision_mask.sum()

		if n_collisions > 0:
			# Remove polymerases that are projected to collide with replisomes
			activeRnaPolys.delByIndexes(np.where(collision_mask)[0])

			# Increment counts of inactive RNA polymerases
			self.inactiveRnaPolys.countInc(n_collisions)

			# Get lengths of transcripts that were terminated prematurely as a
			# result of the collision
			incomplete_sequence_lengths = updated_lengths[collision_mask]

			# Increment counts of bases in incomplete transcripts
			incomplete_sequences = buildSequences(
				self.rnaSequences, TU_indexes[collision_mask],
				np.zeros(n_collisions, dtype=np.int64),
				incomplete_sequence_lengths.max())

			base_counts = np.zeros(4, dtype=np.int64)

			for sl, seq in izip(incomplete_sequence_lengths, incomplete_sequences):
				base_counts += np.bincount(seq[:sl], minlength = 4)

			self.fragmentBases.countsInc(base_counts)
			self.ppi.countInc(n_collisions)

		# Determine if transcript has reached the end of the sequence
		terminalLengths = self.rnaLengths[TU_indexes]
		didTerminate = np.logical_and(
			updated_lengths == terminalLengths,	~collision_mask)
		terminatedRnas = np.bincount(
			TU_indexes[didTerminate], minlength = self.rnaSequences.shape[0])

		# Remove polymerases that have finished transcription from unique
		# molecules
		activeRnaPolys.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = ntpsUsed.sum()

		# Update bulk molecule counts
		self.ntps.countsDec(ntpsUsed)
		self.bulkRnas.countsIs(terminatedRnas)
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

		self.writeToListener("RnapData", "n_collisions", n_collisions)
