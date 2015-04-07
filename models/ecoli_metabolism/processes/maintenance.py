#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

class Maintenance(wholecell.processes.process.Process):
	""" Maintenance """

	_name = "Maintenance"

	# Constructor
	def __init__(self):
		super(Maintenance, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Maintenance, self).initialize(sim, kb)

		# Load constants

		self.cellCycleLen = kb.doubling_time.asNumber(units.s)

		self.initialMaintenanceReactions = (
			kb.NGAM * kb.constants.nAvogadro * kb.mass.avgCellDryMassInit
			).asNumber(1 / units.s) * self.timeStepSec

		# Create views on state
		self.reactants = self.bulkMoleculesView(["ATP[c]", "H2O[c]"])
		self.products = self.bulkMoleculesView(["ADP[c]", "PI[c]", "H[c]"])


	def calculateRequest(self):
		nReactions = self.initialMaintenanceReactions * np.exp(
			np.log(2) / self.cellCycleLen * self.time()
			)

		self.reactants.requestIs(nReactions)


	def evolveState(self):
		reactantCounts = self.reactants.counts()

		assert np.all(reactantCounts == reactantCounts[0])

		self.products.countsInc(reactantCounts[0])

		self.reactants.countsIs(0)
