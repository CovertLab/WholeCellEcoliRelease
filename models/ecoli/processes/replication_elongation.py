#!/usr/bin/env python

"""
Replication

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import polymerize, PAD_VALUE
from wholecell.utils import units

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		# Load parameters
		self.dnaPolymeraseElongationRate = kb.constants.dnaPolymeraseElongationRate.asNumber(units.nt / units.s) * self.timeStepSec
		self.dnaPolymeraseElongationRate = int(round(self.dnaPolymeraseElongationRate)) # TODO: Make this not a hack in the KB

		self.sequenceLengths = kb.process.replication.sequence_lengths
		self.sequences = kb.process.replication.replicaiton_sequences
		self.ntWeights = kb.process.transcription.replicationMonomerWeights

		# Views
		self.activeDnaPoly = self.uniqueMoleculesView('dnaPolymerase')

		self.dntps = self.bulkMoleculesView(kb.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

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
		dNtpCounts = self.dntps.counts()

		activeDnaPoly = self.activeDnaPoly.molecules()

		if len(activeDnaPoly) == 0:
			return

		sequenceIdx, sequeneLengths, massDiffDna = activeDnaPoly.attrs(
			'sequenceIdx', 'sequenceLength', 'massDiff_DNA'
			)

		dNtpsUsed = np.zeros_like(dNtpCounts)

		sequences = buildSequences(
			self.sequences,
			sequenceIdx,
			sequeneLengths,
			self.dnaPolymeraseElongationRate
			)

		reactionLimit = dNtpCounts.sum() # TODO: account for energy

		sequenceElongations, dNtpsUsed, nElongations = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		massIncreaseDna = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.ntWeights.asUnit(units.fg)
			)

		updatedMass = massDiffDna + massIncreaseDna

		didInitialize = (sequeneLengths == 0) & (sequenceElongations > 0)

		updatedLengths = sequeneLengths + sequenceElongations

		activeDnaPoly.attrIs(
			sequenceLength = updatedLengths,
			massDiff_DNA = updatedMass
			)

		terminalLengths = self.sequence_lengths[sequenceIdx]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedChromosomes = np.bincount(
			sequenceIdx[didTerminate],
			minlength = self.sequences.shape[0]
			)

		activeDnaPoly.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = dNtpsUsed.sum()

		self.dntps.countsDec(dNtpsUsed)

		self.bulkChromosomes.countsIs(terminatedChromosomes)

		self.inactiveRnaPolys.countInc(nTerminated)

		self.ppi.countInc(nElongations - nInitialized)

		# expectedElongations = np.fmin(
		# 	self.elngRate,
		# 	terminalLengths - sequeneLengths
		# 	)

		#dnaPolymeraseStalls = expectedElongations - sequenceElongations
