#!/usr/bin/env python

"""
ReplicationElongation

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease, PAD_VALUE
from wholecell.utils import units

class ReplicationElongation(wholecell.processes.process.Process):
	""" ReplicationElongation """

	_name = "ReplicationElongation"

	# Constructor
	def __init__(self):
		super(ReplicationElongation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ReplicationElongation, self).initialize(sim, sim_data)

		# Load parameters
		self.criticalInitiationMass = sim_data.growthRateParameters.getDnaCriticalMass(sim_data.conditionToDoublingTime[sim_data.condition])
		self.getDnaCriticalMass = sim_data.growthRateParameters.getDnaCriticalMass
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		self.sequenceLengths = sim_data.process.replication.sequence_lengths
		self.sequences = sim_data.process.replication.replication_sequences
		self.polymerized_dntp_weights = sim_data.process.replication.replicationMonomerWeights

		self.dnaPolyElngRate = int(round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(units.nt / units.s)))

		# Views
		self.activeDnaPoly = self.uniqueMoleculesView('dnaPolymerase')

		self.oriCs = self.uniqueMoleculesView('originOfReplication')

		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.chromosomeHalves = self.bulkMoleculesView(sim_data.moleculeGroups.partialChromosome)

		self.full_chromosome = self.bulkMoleculeView("CHROM_FULL[c]")

	def calculateRequest(self):
		self.criticalInitiationMass = self.getDnaCriticalMass(self.nutrientToDoublingTime[self._sim.processes["PolypeptideElongation"].currentNutrients])

		self.full_chromosome.requestAll()
		self.oriCs.requestAll()

		activeDnaPoly = self.activeDnaPoly.allMolecules()

		if len(activeDnaPoly) == 0:
			return

		self.activeDnaPoly.requestAll()
		
		sequenceIdx, sequenceLength = activeDnaPoly.attrs(
			'sequenceIdx', 'sequenceLength'
			)

		sequences = buildSequences(
			self.sequences,
			sequenceIdx,
			sequenceLength,
			self._dnaPolymeraseElongationRate()
			)

		sequenceComposition = np.bincount(sequences[sequences != PAD_VALUE], minlength = 4)

		dNtpsTotal = self.dntps.total()

		maxFractionalReactionLimit = (np.fmin(1, dNtpsTotal/sequenceComposition)).min()

		self.dntps.requestIs(
			maxFractionalReactionLimit * sequenceComposition
			)

	# Calculate temporal evolution
	def evolveState(self):

		activeDnaPoly = self.activeDnaPoly.molecules()

		##########################################
		# Perform replication initiation process #
		##########################################

		activeDnaPoly = self.activeDnaPoly.molecules()
		activePolymerasePresent = len(activeDnaPoly) > 0

		oriCs = self.oriCs.molecules()
		chromosomes = self.full_chromosome.total()[0]

		if len(oriCs) == 0 and chromosomes == 0:
			return

		if activePolymerasePresent:
			replicationRound = activeDnaPoly.attr('replicationRound')

		# Get cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Initiate if over a critical mass threshold, and no oriC currently exists which was initiated at that same critical mass.
		initiate = False
		massFactor = cellMass / self.criticalInitiationMass
		massPerOrigin = massFactor / len(oriCs)

		self.writeToListener("ReplicationData", "criticalMassPerOriC", massPerOrigin)
		self.writeToListener("ReplicationData", "criticalInitiationMass", self.criticalInitiationMass.asNumber(units.fg))

		if massPerOrigin >= 1.0:
			initiate = True

		if initiate:
			# Number of oriC the cell has
			if activePolymerasePresent:
				numOric = 2 * np.unique(replicationRound).size * self.full_chromosome.count()
			else:
				numOric = 1 * self.full_chromosome.count()
				replicationRound = np.array([0])
		
			numberOfNewPolymerase = 4 * numOric

			# Initialize 4 DNA polymerases per oriC
			activeDnaPoly = self.activeDnaPoly.moleculesNew(
				"dnaPolymerase",
				numberOfNewPolymerase
				)

			# Generate a new oriC for each old oriC
			oriCs = self.oriCs.moleculesNew(
				"originOfReplication",
				numOric
				)

			sequenceIdx = np.tile(np.array([0,1,2,3], dtype=np.int8), numOric)
			sequenceLength = np.zeros(numberOfNewPolymerase, dtype = np.int8)
			replicationRound = np.ones(numberOfNewPolymerase, dtype=np.int8) * (replicationRound.max() + 1)
			chromosomeIndex = np.zeros(numberOfNewPolymerase, dtype=np.int8)
			chromosomeIndex[numberOfNewPolymerase / 2:] = 1.

			activeDnaPoly.attrIs(
				sequenceIdx = sequenceIdx,
				sequenceLength = sequenceLength,
				replicationRound = replicationRound,
				chromosomeIndex = chromosomeIndex,
				)

		##########################################
		# Perform replication elongation process #
		##########################################

		if len(activeDnaPoly) == 0:
			return
		
		dNtpCounts = self.dntps.counts()

		sequenceIdx, sequenceLengths, massDiffDna = activeDnaPoly.attrs(
			'sequenceIdx', 'sequenceLength', 'massDiff_DNA'
			)

		dNtpsUsed = np.zeros_like(dNtpCounts)

		sequences = buildSequences(
			self.sequences,
			sequenceIdx,
			sequenceLengths,
			self._dnaPolymeraseElongationRate()
			)

		reactionLimit = dNtpCounts.sum() # TODO: account for energy

		sequenceElongations, dNtpsUsed, nElongations = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		# Compute mass increase for each polymerizing chromosome half
		massIncreaseDna = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.polymerized_dntp_weights.asNumber(units.fg)
			)

		updatedMass = massDiffDna + massIncreaseDna

		didInitialize = (sequenceLengths == 0) & (sequenceElongations > 0)

		updatedLengths = sequenceLengths + sequenceElongations

		activeDnaPoly.attrIs(
			sequenceLength = updatedLengths,
			massDiff_DNA = updatedMass
			)

		# Determine if any chromosome half finshed polymerizing
		# and create new unique chromosomes
		terminalLengths = self.sequenceLengths[sequenceIdx]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedChromosomes = np.bincount(
			sequenceIdx[didTerminate],
			minlength = self.sequences.shape[0]
			)

		if didTerminate.sum():
			# If any forks terminated, update the chromosome index of the next fork up the chromosome
			sequenceIdx, sequenceLength, chromosomeIndex, replicationRound = activeDnaPoly.attrs(
				'sequenceIdx', 'sequenceLength', 'chromosomeIndex', 'replicationRound'
				)

			# Serach for next fork up the chromosome
			# Criteria: Must match the same sequence index, must be current terminating forks's replication
			# round +1, must match terminating fork's chromosome index (so you don't flip the same upstream
			# fork twice).
			sequenceIdxMatch = np.zeros(didTerminate.shape)
			replicationRoundMatch = np.zeros(didTerminate.shape)
			chromosomeIndexMatch = np.zeros(didTerminate.shape)
			for idx in np.where(didTerminate)[0]:
				sequenceIdxMatch = np.logical_or(sequenceIdxMatch, sequenceIdx[idx] == sequenceIdx)
				replicationRoundMatch = np.logical_or(replicationRoundMatch, replicationRound[idx] + 1 == replicationRound)
				chromosomeIndexMatch = np.logical_or(chromosomeIndexMatch, chromosomeIndex[idx] == chromosomeIndex)

			# Potential matches for chromosome index to flip
			potentialMatches = np.logical_and.reduce((sequenceIdxMatch, replicationRoundMatch, chromosomeIndexMatch))

			# Take the first two. If there are more they are degenerate.
			toIndex = np.where(potentialMatches)[0][:2]

			# Changes 0 -> 1 and 1 -> 0
			toFlip = chromosomeIndex[toIndex]
			chromosomeIndex[toIndex] = np.abs(toFlip - 1)

			# Resets chromosomeIndex
			activeDnaPoly.attrIs(chromosomeIndex = chromosomeIndex)

		newUniqueChromosomeHalves = sequenceIdx[np.where(didTerminate)[0]]

		activeDnaPoly.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = dNtpsUsed.sum()

		self.chromosomeHalves.countsInc(terminatedChromosomes)

		self.dntps.countsDec(dNtpsUsed)

		self.ppi.countInc(nElongations)

	def _dnaPolymeraseElongationRate(self):
		return self.dnaPolyElngRate * self.timeStepSec()
