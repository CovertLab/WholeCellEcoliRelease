#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound

class TranscriptElongation(wholecell.processes.process.Process):
	""" TranscriptElongation """

	_name = "TranscriptElongation"

	def __init__(self):
		super(TranscriptElongation, self).__init__()


	def initialize(self, sim, kb):
		super(TranscriptElongation, self).initialize(sim, kb)

		# TODO: account for the polymer ends

		# Load constants

		## Molecule IDs

		ntpIDs = kb.ntpIds

		## Find the average RNA composition

		synthProb = kb.rnaData["synthProb"].magnitude
		compositionAll = kb.rnaData["countsACGU"].magnitude

		# TODO: better model the variance of this distribution

		compositionUnnormed = np.dot(compositionAll.T, synthProb)

		self.monomerComposition = compositionUnnormed / compositionUnnormed.sum()

		## Find the average total transcription rate with respect to cell age

		self.cellCycleLen = kb.cellCycleLen.to("s").magnitude

		initialDryMass = kb.avgCellDryMassInit.to("fg").magnitude

		rnaMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].to("min").magnitude == 60.0
			]["rnaMassFraction"]

		initialRnaMass = initialDryMass * rnaMassFraction

		monomerMWs = kb.transcriptionMonomerWeights

		monomerAverageMW = np.dot(monomerMWs, self.monomerComposition) # average MW weighted by transcript composition

		self.initialAverageMonomerCounts = initialRnaMass / monomerAverageMW

		# Create views on state

		self.ntps = self.bulkMoleculesView(ntpIDs)
		self.ppi = self.bulkMoleculeView("PPI[c]")

		# TODO: incorporated nucleotides


	def calculateRequest(self):
		totalMonomers = np.int64(stochasticRound(
			self.randomState,
			self.initialAverageMonomerCounts
			* np.exp(np.log(2) / self.cellCycleLen * self.time())
			* (np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1)
			))

		ntpsRequested = self.randomState.multinomial(
			totalMonomers,
			self.monomerComposition
			)

		if (ntpsRequested > self.ntps.total()).any():
			# TODO: flag simulation instead of printing
			print "{} is metabolically limited".format(self.name())

		self.ntps.requestIs(ntpsRequested)


	def evolveState(self):
		self.ppi.countInc(self.ntps.counts().sum())

		self.ntps.countsIs(0)
