#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# TODO: account for the polymer ends

		# Load constants

		## Molecule IDs

		aaIDs = kb.aaIDs

		## Find the average RNA composition

		synthProbUnnormed = kb.rnaExpression["expression"].magnitude[kb.rnaIndexToMonomerMapping]

		synthProb = synthProbUnnormed / synthProbUnnormed.sum()
		compositionAll = kb.monomerData["aaCounts"].magnitude

		# TODO: better model the variance of this distribution

		compositionUnnormed = np.dot(compositionAll.T, synthProb)

		self.monomerComposition = compositionUnnormed / compositionUnnormed.sum()

		## Find the average total transcription rate with respect to cell age

		self.cellCycleLen = kb.cellCycleLen.to("s").magnitude

		initialDryMass = kb.avgCellDryMassInit.to("fg").magnitude

		proteinMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].to("min").magnitude == 60.0
			]["proteinMassFraction"]

		initialProteinMass = initialDryMass * proteinMassFraction

		monomerMWs = kb.translationMonomerWeights

		monomerAverageMW = np.dot(monomerMWs, self.monomerComposition) # average MW weighted by transcript composition

		self.initialAverageMonomerCounts = initialProteinMass / monomerAverageMW

		## Energy costs

		self.gtpPerElongation = kb.gtpPerTranslation

		# Create views on state

		self.aas = self.bulkMoleculesView(aaIDs)
		self.h2o = self.bulkMoleculeView("H2O[c]")
		self.gtp = self.bulkMoleculeView("GTP[c]")

		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h = self.bulkMoleculeView("H[c]")

		# TODO: incorporated AAs


	def calculateRequest(self):
		totalMonomers = np.int64(stochasticRound(
			self.randomState,
			self.initialAverageMonomerCounts
			* np.exp(np.log(2) / self.cellCycleLen * self.time())
			* (np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1)
			))

		aasRequested = self.randomState.multinomial(
			totalMonomers,
			self.monomerComposition
			)

		gtpRequested = np.int64(stochasticRound(
			self.randomState,
			self.gtpPerElongation * aasRequested.sum()
			))

		self.aas.requestIs(aasRequested)
		self.gtp.requestIs(gtpRequested)
		self.h2o.requestIs(gtpRequested) # NOTE: this is an overestimate


	def evolveState(self):
		self.h2o.countInc(self.aas.counts().sum())

		self.gdp.countInc(self.gtp.count())
		self.pi.countInc(self.gtp.count())
		self.h.countInc(self.gtp.count())

		self.h2o.countDec(self.gtp.count())

		self.aas.countsIs(0)
		self.gtp.countIs(0)
