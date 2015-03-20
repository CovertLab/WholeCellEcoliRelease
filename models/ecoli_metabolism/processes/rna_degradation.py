#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

from wholecell.utils.fitting import countsFromMassAndExpression

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		# Load constants

		## Molecule IDs

		nmpIDs = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"]
		polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedNT_IDs]

		## Find the magnitude and distribution of amino acids recovered by degradation

		self.cellCycleLen = kb.constants.cellCycleLen.asNumber(units.s)

		rnaComposition = kb.process.transcription.rnaData["countsACGU"].asNumber()

		initialDryMass = kb.constants.avgCellDryMassInit

		rnaMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].asNumber(units.min) == 60.0
			]["rnaMassFraction"]

		initialRnaMass = initialDryMass * rnaMassFraction

		initialRnaCounts = countsFromMassAndExpression(
			initialRnaMass.asNumber(units.g),
			kb.process.transcription.rnaData["mw"].asNumber(units.g / units.mol),
			kb.process.transcription.rnaData["expression"],
			kb.constants.nAvogadro.asNumber(1 / units.mol)
			) * kb.process.transcription.rnaData["expression"]

		initialRnaDegradationRate = initialRnaCounts * (
			kb.process.transcription.rnaData["degRate"]
			).asNumber(1 / units.s)

		initialDegrading = np.dot(rnaComposition.T, initialRnaDegradationRate)

		self.initialDegradingTotal = initialDegrading.sum()

		self.monomerComposition = initialDegrading / initialDegrading.sum()

		# Create views on state

		self.polymerized = self.bulkMoleculesView(polymerizedIDs)
		self.h2o = self.bulkMoleculeView("H2O[c]")

		self.nmps = self.bulkMoleculesView(nmpIDs)
		self.proton = self.bulkMoleculeView("H[c]")


	def calculateRequest(self):
		totalMonomers = np.int64(stochasticRound(
			self.randomState,
			self.initialDegradingTotal
			* np.exp(np.log(2) / self.cellCycleLen * self.time())
			))

		polymerizedRequested = self.randomState.multinomial(
			totalMonomers,
			self.monomerComposition
			)

		h2oRequested = polymerizedRequested.sum()

		if h2oRequested > self.h2o.total():

			# TODO: flag simulation instead of printing
			print "{} is metabolically limited".format(self.name())

		self.polymerized.requestIs(polymerizedRequested)

		self.h2o.requestIs(h2oRequested)


	def evolveState(self):
		polymerizedCounts = self.polymerized.counts()

		self.nmps.countsInc(polymerizedCounts)
		self.proton.countInc(polymerizedCounts.sum())

		self.h2o.countDec(polymerizedCounts.sum())

		self.polymerized.countsIs(0)
