#!/usr/bin/env python

"""
UniqueTranscriptElongation

Transcription elongation sub-model.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize_new import buildSequences, polymerize, computeMassIncrease, PAD_VALUE

# TODO: refactor mass calculations
# TODO: confirm reaction stoich
# TODO: resolve mounting process namespace issues

class UniqueTranscriptElongation(wholecell.processes.process.Process):
	""" UniqueTranscriptElongation """

	_name = "UniqueTranscriptElongation"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None
		self.rnaIds = None
		self.rnaLengths = None
		self.rnaSequences = None
		self.ntWeights = None
		self.hydroxylWeight = None

		# Views
		self.activeRnaPolys = None
		self.bulkRnas = None
		self.ntps = None
		self.ppi = None
		self.h2o = None
		self.proton = None
		self.rnapSubunits = None

		# Cached values
		self._sequences = None
		self._transcriptLengths = None
		self._rnaIndexes = None

		super(UniqueTranscriptElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniqueTranscriptElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = kb.rnaPolymeraseElongationRate.to('nucleotide / s').magnitude

		self.rnaIds = kb.rnaData['id']

		# TODO: refactor mass updates

		sequences = kb.rnaData["sequence"]

		self.rnaLengths = kb.rnaData["length"].magnitude

		maxLen = np.int64(self.rnaLengths.max() + self.elngRate)

		self.rnaSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.rnaSequences.fill(PAD_VALUE)

		ntMapping = {ntpId:i for i, ntpId in enumerate(["A", "C", "G", "U"])}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.rnaSequences[i, j] = ntMapping[letter]

		# TOKB
		self.ntWeights = np.array([
			345.20, # A
			321.18, # C
			361.20, # G
			322.17, # U
			]) - 17.01 # weight of a hydroxyl

		# TOKB
		self.hydroxylWeight = 17.01 # counted once for the end of the polymer

		self.ntWeights *= 1e15/6.022e23
		self.hydroxylWeight *= 1e15/6.022e23

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)

		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")


	def calculateRequest(self):
		activeRnaPolys = self.activeRnaPolys.allMolecules()

		if len(activeRnaPolys) == 0:
			return

		self.activeRnaPolys.requestAll()

		rnaIndexes, transcriptLengths = activeRnaPolys.attrs(
			'rnaIndex', 'transcriptLength'
			)

		sequences = buildSequences(
			self.rnaSequences,
			rnaIndexes,
			transcriptLengths,
			self.elngRate
			)

		sequenceComposition = np.bincount(sequences[sequences != PAD_VALUE], minlength = 4)

		ntpsTotal = self.ntps.total()

		maxFractionalReactionLimit = (np.fmin(1, ntpsTotal/sequenceComposition)).min()

		self.ntps.requestIs(
			maxFractionalReactionLimit * sequenceComposition
			)

		self.h2o.requestIs(self.ntps.total().sum()) # this drastically overestimates water assignment


	# Calculate temporal evolution
	def evolveState(self):
		ntpCounts = self.ntps.counts()

		activeRnaPolys = self.activeRnaPolys.molecules()

		if len(activeRnaPolys) == 0:
			return

		rnaIndexes, transcriptLengths, massDiffRna = activeRnaPolys.attrs(
			'rnaIndex', 'transcriptLength', 'massDiffRna'
			)

		ntpsUsed = np.zeros_like(ntpCounts)

		sequences = buildSequences(
			self.rnaSequences,
			rnaIndexes,
			transcriptLengths,
			self.elngRate
			)

		reactionLimit = ntpCounts.sum() # TODO: account for energy

		sequenceElongations, ntpsUsed, nElongations = polymerize(
			sequences,
			ntpCounts,
			reactionLimit,
			self.randomState
			)

		massIncreaseRna = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.ntWeights
			)

		updatedMass = massDiffRna + massIncreaseRna

		didInitialize = (transcriptLengths == 0) & (sequenceElongations > 0)

		updatedLengths = transcriptLengths + sequenceElongations

		updatedMass[didInitialize] += self.hydroxylWeight

		activeRnaPolys.attrIs(
			transcriptLength = updatedLengths,
			massDiffRna = updatedMass
			)

		didTerminate = (updatedLengths == self.rnaLengths[rnaIndexes])

		terminatedRnas = np.bincount(
			rnaIndexes[didTerminate],
			minlength = self.rnaSequences.shape[0]
			)

		activeRnaPolys.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = ntpsUsed.sum()

		self.ntps.countsDec(ntpsUsed)

		self.bulkRnas.countsIs(terminatedRnas)

		self.inactiveRnaPolys.countInc(nTerminated)

		self.h2o.countDec(nInitialized)
		self.proton.countInc(nInitialized)

		self.ppi.countInc(nElongations)
