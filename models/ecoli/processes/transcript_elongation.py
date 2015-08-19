#!/usr/bin/env python

"""
TranscriptElongation

Transcription elongation sub-model.

TODO:
- use transcription units instead of single genes

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease, PAD_VALUE
from wholecell.utils import units

class TranscriptElongation(wholecell.processes.process.Process):
	""" TranscriptElongation """

	_name = "TranscriptElongation"

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

		super(TranscriptElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(TranscriptElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = kb.growthRateParameters.rnaPolymeraseElongationRate.asNumber(units.nt / units.s) * self.timeStepSec
		self.elngRate = int(round(self.elngRate)) # TODO: Make this less of a hack by implementing in the KB

		self.rnaIds = kb.process.transcription.rnaData['id']

		self.rnaLengths = kb.process.transcription.rnaData["length"].asNumber()

		self.rnaSequences = kb.process.transcription.transcriptionSequences

		self.ntWeights = kb.process.transcription.transcriptionMonomerWeights

		self.endWeight = kb.process.transcription.transcriptionEndWeight

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)

		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')

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

		self.writeToListener("GrowthLimits", "ntpPoolSize", self.ntps.total())
		self.writeToListener("GrowthLimits", "ntpRequestSize", maxFractionalReactionLimit * sequenceComposition)

	# Calculate temporal evolution
	def evolveState(self):
		ntpCounts = self.ntps.counts()

		self.writeToListener("GrowthLimits", "ntpAllocated", self.ntps.counts())

		activeRnaPolys = self.activeRnaPolys.molecules()

		if len(activeRnaPolys) == 0:
			return

		rnaIndexes, transcriptLengths, massDiffRna = activeRnaPolys.attrs(
			'rnaIndex', 'transcriptLength', 'massDiff_mRNA'
			)

		ntpsUsed = np.zeros_like(ntpCounts)

		sequences = buildSequences(
			self.rnaSequences,
			rnaIndexes,
			transcriptLengths,
			self.elngRate
			)

		ntpCountInSequence = np.bincount(sequences[sequences != PAD_VALUE], minlength = 4)

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

		updatedMass[didInitialize] += self.endWeight

		activeRnaPolys.attrIs(
			transcriptLength = updatedLengths,
			massDiff_mRNA = updatedMass
			)

		terminalLengths = self.rnaLengths[rnaIndexes]

		didTerminate = (updatedLengths == terminalLengths)

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

		self.ppi.countInc(nElongations - nInitialized)

		expectedElongations = np.fmin(
			self.elngRate,
			terminalLengths - transcriptLengths
			)

		rnapStalls = expectedElongations - sequenceElongations

		self.writeToListener("GrowthLimits", "ntpUsed", ntpsUsed)

		self.writeToListener("RnapData", "rnapStalls", rnapStalls)
		self.writeToListener("RnapData", "ntpCountInSequence", ntpCountInSequence)
		self.writeToListener("RnapData", "ntpCounts", ntpCounts)

		self.writeToListener("RnapData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RnapData", "actualElongations", sequenceElongations.sum())

		self.writeToListener("RnapData", "didTerminate", didTerminate.sum())
		self.writeToListener("RnapData", "terminationLoss", (terminalLengths - transcriptLengths)[didTerminate].sum())
