#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

from wholecell.utils.fitting import calcProteinCounts

class ProteinDegradation(wholecell.processes.process.Process):
	""" ProteinDegradation """

	_name = "ProteinDegradation"

	# Constructor
	def __init__(self):
		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ProteinDegradation, self).initialize(sim, kb)

		# Load constants

		## Molecule IDs

		aaIDs = kb.moleculeGroups.aaIDs
		polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedAA_IDs]

		## Find the magnitude and distribution of amino acids recovered by degradation

		self.cellCycleLen = kb.constants.cellCycleLen.asNumber(units.s)

		proteinComposition = kb.process.translation.monomerData["aaCounts"].asNumber()

		initialDryMass = kb.constants.avgCellDryMassInit

		proteinMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].asNumber(units.min) == 60.0
			]["proteinMassFraction"]

		initialProteinMass = initialDryMass * proteinMassFraction

		initialProteinCounts = calcProteinCounts(kb, initialProteinMass)

		initialProteinDegradationRate = initialProteinCounts * (
			kb.process.translation.monomerData["degRate"] # NOTE: constrast this with the translation submodel, which accounts for degradation AND dilution
			).asNumber(1 / units.s)

		initialDegrading = np.dot(proteinComposition.T, initialProteinDegradationRate)

		self.initialDegradingTotal = initialDegrading.sum()

		self.monomerComposition = initialDegrading / initialDegrading.sum()

		# Create views on state

		self.polymerized = self.bulkMoleculesView(polymerizedIDs)
		self.h2o = self.bulkMoleculeView("H2O[c]")

		self.aas = self.bulkMoleculesView(aaIDs)


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

		self.aas.countsInc(polymerizedCounts)

		self.h2o.countDec(polymerizedCounts.sum())

		self.polymerized.countsIs(0)
