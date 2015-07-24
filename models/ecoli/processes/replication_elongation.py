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
	def initialize(self, sim, kb):
		super(ReplicationElongation, self).initialize(sim, kb)

		# Load parameters
		self.dnaPolymeraseElongationRate = kb.constants.dnaPolymeraseElongationRate.asNumber(units.nt / units.s) * self.timeStepSec
		self.dnaPolymeraseElongationRate = int(round(self.dnaPolymeraseElongationRate)) # TODO: Make this not a hack in the KB

		self.sequenceLengths = kb.process.replication.sequence_lengths
		self.sequences = kb.process.replication.replication_sequences
		self.polymerized_dntp_weights = kb.process.replication.replicationMonomerWeights

		# Views
		self.activeDnaPoly = self.uniqueMoleculesView('dnaPolymerase')

		self.dntps = self.bulkMoleculesView(kb.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.chromosomeHalves = self.bulkMoleculesView(kb.moleculeGroups.partialChromosome)

	def calculateRequest(self):
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
			self.dnaPolymeraseElongationRate
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

		# TODO: Figure out what to do if the chromosome has already divided but more initialization is required!

		activeDnaPoly = self.activeDnaPoly.molecules()
		activePolymerasePresent = len(activeDnaPoly) > 0

		if activePolymerasePresent:
			sequenceLength, replicationRound = activeDnaPoly.attrs('sequenceLength', 'replicationRound')

		# Get cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		refractionOver = True
		if activePolymerasePresent:
			refractionOver = sequenceLength.min() > 5000

		initiate = False
		diffFactor = 0.2
		# TODO: Might need to vary diff factor by growth rate
		if (cellMass > self.criticalInitiationMass) and np.abs((cellMass - self.criticalInitiationMass).asNumber(units.fg)) < diffFactor and refractionOver:
			initiate = True

		if initiate:
			# Number of oriC the cell has
			if activePolymerasePresent:
				numOric = 2 * np.unique(replicationRound).size
			else:
				numOric = 1
				replicationRound = np.array([0])
			
			numberOfNewPolymerase = 4 * numOric

			activeDnaPoly = self.activeDnaPoly.moleculesNew(
				"dnaPolymerase",
				numberOfNewPolymerase
				)

			sequenceIdx = np.tile(np.array([0,1,2,3], dtype=np.int8), numOric)
			sequenceLength = np.zeros(numberOfNewPolymerase, dtype = np.int8)
			replicationRound = np.ones(numberOfNewPolymerase, dtype=np.int8) * (replicationRound.max() + 1)
			replicationDivision = np.zeros(numberOfNewPolymerase, dtype=np.int8)
			replicationDivision[numberOfNewPolymerase / 2:] = 1.

			activeDnaPoly.attrIs(
				sequenceIdx = sequenceIdx,
				sequenceLength = sequenceLength,
				replicationRound = replicationRound,
				replicationDivision = replicationDivision,
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
			self.dnaPolymeraseElongationRate
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

		newUniqueChromosomeHalves = sequenceIdx[np.where(didTerminate)[0]]

		activeDnaPoly.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = dNtpsUsed.sum()

		self.chromosomeHalves.countsInc(terminatedChromosomes)

		self.dntps.countsDec(dNtpsUsed)

		self.ppi.countInc(nElongations - nInitialized)