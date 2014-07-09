#!/usr/bin/env python

"""
AtpUsage

Hydrolyze ATP (for growth associated maintenance)

TODO:
- scale with a macro cell property (mass or volume) instead of time

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_ATP_USAGE

class AtpUsage(wholecell.processes.process.Process):
	""" AtpUsage """

	_name = "AtpUsage"

	# Constructor
	def __init__(self):
		super(AtpUsage, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(AtpUsage, self).initialize(sim, kb)

		# Load constants
		self.nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
		self.initialDryMass = kb.avgCellDryMassInit.to('g').magnitude
		self.cellCycleLen = kb.cellCycleLen.to('s').magnitude

		moleculeIds = ["ATP[c]", "H2O[c]", "PI[c]", "ADP[c]", "H[c]"]
		self.molecules = self.bulkMoleculesView(moleculeIds)
		self.reactants = self.bulkMoleculesView(["ATP[c]", "H2O[c]"])
		self.products = self.bulkMoleculesView(["ADP[c]", "PI[c]", "H[c]"])
		self.atp = self.bulkMoleculeView("ATP[c]")
		self.h2o = self.bulkMoleculeView("H2O[c]")
		# self.pi = self.bulkMoleculeView("PI[c]")
		# self.adp = self.bulkMoleculeView("ADP[c]")
		# self.h = self.bulkMoleculeView("H[c]")

		self.atpInitialPool = (
			kb.atpPoolSize * self.timeStepSec *
			kb.nAvogadro.to("1/millimole") * kb.avgCellDryMassInit.to("DCW_gram")
			).magnitude

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_ATP_USAGE)

	def calculateRequest(self):

		poolToRequest = self.atpInitialPool * np.exp(np.log(2) / self.cellCycleLen * self.time())
		self.reactants.requestIs(poolToRequest)

	def evolveState(self):
		
		atpsHydrolyzed = np.fmin(
			self.atp.count(),
			self.h2o.count()
			)
		self.reactants.countsDec(atpsHydrolyzed)
		self.products.countsInc(atpsHydrolyzed)
