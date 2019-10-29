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
		self.transcription_data = sim_data.process.transcription

		# Views
		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)
		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.flat_elongation = not sim._variable_elongation_transcription

	def calculateRequest(self):
		# Calculate elongation rate based on the current nutrients
		current_nutrients = self._external_states['Environment'].nutrients

		self.rnapElongationRate = self.rnaPolymeraseElongationRateDict[current_nutrients].asNumber(units.nt / units.s)

		self.elongation_rates = self.transcription_data.make_elongation_rates(
			self.randomState,
			self.rnapElongationRate,
			self.timeStepSec(),
			self.flat_elongation)

		# Request all active RNA polymerases
		activeRnaPolys = self.activeRnaPolys.allMolecules()
		if len(activeRnaPolys) == 0:
			return
		self.activeRnaPolys.requestAll()

		# Determine total possible sequences of nucleotides that can be transcribed in this time step for each polymerase
		rnaIndexes, transcriptLengths = activeRnaPolys.attrs('rnaIndex', 'transcriptLength')
		sequences = buildSequences(
			self.rnaSequences,
			rnaIndexes,
			transcriptLengths,
			self.elongation_rates)

		sequenceComposition = np.bincount(sequences[sequences != polymerize.PAD_VALUE], minlength = 4)

		# Calculate if any nucleotides are limited and request up to the number in the sequences or number available
		ntpsTotal = self.ntps.total()
		maxFractionalReactionLimit = np.fmin(1, ntpsTotal / sequenceComposition)
		self.ntps.requestIs(maxFractionalReactionLimit * sequenceComposition)

		self.writeToListener("GrowthLimits", "ntpPoolSize", self.ntps.total())
		self.writeToListener("GrowthLimits", "ntpRequestSize", maxFractionalReactionLimit * sequenceComposition)

	def evolveState(self):
		ntpCounts = self.ntps.counts()
		self.writeToListener("GrowthLimits", "ntpAllocated", ntpCounts)

		activeRnaPolys = self.activeRnaPolys.molecules()
		if len(activeRnaPolys) == 0:
			return

		# Determine sequences that can be elongated
		rnaIndexes, transcriptLengths, massDiffRna = activeRnaPolys.attrs('rnaIndex', 'transcriptLength', 'massDiff_mRNA')
		sequences = buildSequences(
			self.rnaSequences,
			rnaIndexes,
			transcriptLengths,
			self.elongation_rates)

		ntpCountInSequence = np.bincount(sequences[sequences != polymerize.PAD_VALUE], minlength = 4)

		# Polymerize transcripts based on sequences and available nucleotides
		reactionLimit = ntpCounts.sum()
		active_elongation_rates = self.elongation_rates[rnaIndexes]

		result = polymerize(sequences, ntpCounts, reactionLimit, self.randomState, active_elongation_rates)
		sequenceElongations = result.sequenceElongation
		ntpsUsed = result.monomerUsages
		nElongations = result.nReactions

		# Calculate changes in mass associated with polymerization and update active polymerases
		massIncreaseRna = computeMassIncrease(sequences, sequenceElongations, self.ntWeights)
		updatedMass = massDiffRna + massIncreaseRna
		didInitialize = (transcriptLengths == 0) & (sequenceElongations > 0)
		updatedLengths = transcriptLengths + sequenceElongations
		updatedMass[didInitialize] += self.endWeight
		activeRnaPolys.attrIs(transcriptLength = updatedLengths, massDiff_mRNA = updatedMass)

		# Determine if transcript has reached the end of the sequence
		terminalLengths = self.rnaLengths[rnaIndexes]
		didTerminate = (updatedLengths == terminalLengths)
		terminatedRnas = np.bincount(rnaIndexes[didTerminate], minlength = self.rnaSequences.shape[0])

		# Remove polymerases that have finished transcription from unique molecules
		activeRnaPolys.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = ntpsUsed.sum()

		# Update bulk molecule counts
		self.ntps.countsDec(ntpsUsed)
		self.bulkRnas.countsIs(terminatedRnas)
		self.inactiveRnaPolys.countInc(nTerminated)
		self.ppi.countInc(nElongations - nInitialized)

		# Calculate stalls
		expectedElongations = np.fmin([self.elongation_rates[x] for x in rnaIndexes], terminalLengths - transcriptLengths)
		rnapStalls = expectedElongations - sequenceElongations

		# Write outputs to listeners
		self.writeToListener("TranscriptElongationListener", "countRnaSynthesized", terminatedRnas)
		self.writeToListener("TranscriptElongationListener", "countNTPsUSed", nElongations)
		self.writeToListener("GrowthLimits", "ntpUsed", ntpsUsed)
		self.writeToListener("RnapData", "rnapStalls", rnapStalls)
		self.writeToListener("RnapData", "ntpCountInSequence", ntpCountInSequence)
		self.writeToListener("RnapData", "ntpCounts", ntpCounts)
		self.writeToListener("RnapData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RnapData", "actualElongations", sequenceElongations.sum())
		self.writeToListener("RnapData", "didTerminate", didTerminate.sum())
		self.writeToListener("RnapData", "terminationLoss", (terminalLengths - transcriptLengths)[didTerminate].sum())
