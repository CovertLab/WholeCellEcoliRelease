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

import warnings

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_ATP_USAGE
from wholecell.utils.random import stochasticRound

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

		self.growthAssociated_reactionsPerFemtogram = (
			kb.atpUsedPerMassIncrease * kb.nAvogadro
			).to("1/femtogram").magnitude

		self.nongrowthAssociated_reactionsPerTimestep = (
			kb.atpUsedPerSecond * kb.nAvogadro
			).to("1/femtogram/s").magnitude * self.timeStepSec

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_ATP_USAGE)


	def calculateRequest(self):
		mass = self.readFromListener("Mass", "dryMass")
		deltaMass = self.readFromListener("Mass", "growth")

		expectedReactions_growthAssociated = np.fmax(
			deltaMass * self.growthAssociated_reactionsPerFemtogram,
			0
			)

		expectedReactions_nongrowthAssociated = (mass
			* self.nongrowthAssociated_reactionsPerTimestep)

		expectedReactions = (expectedReactions_growthAssociated
			+ expectedReactions_nongrowthAssociated)
		
		self.reactants.requestIs(
			stochasticRound(self.randomState, expectedReactions)
			)


	def evolveState(self):
		mass = self.readFromListener("Mass", "dryMass")
		deltaMass = self.readFromListener("Mass", "growth")

		expectedReactions_growthAssociated = np.fmax(
			deltaMass * self.growthAssociated_reactionsPerFemtogram,
			0
			)

		expectedReactions_nongrowthAssociated = (mass
			* self.nongrowthAssociated_reactionsPerTimestep)

		expectedReactions = (expectedReactions_growthAssociated
			+ expectedReactions_nongrowthAssociated)

		import ipdb; ipdb.set_trace()
		
		atpsHydrolyzed = np.fmin(
			self.atp.count(),
			self.h2o.count()
			)

		self.reactants.countsDec(atpsHydrolyzed)
		self.products.countsInc(atpsHydrolyzed)

		if atpsHydrolyzed < np.floor(expectedReactions):
			warnings.warn("Maintenance not satisfied; the cell may be growing too fast")

			# TODO: flag some sort of listener that tracks phenomenological
			# observations like cell death/division

		# TODO: consider implementing (N)GAM as a metabolism constraint
		# TODO: let (N)GAM roll over to the next time step, just in case
		# there is some stochasticity (i.e. ATP limiting in one time step, in 
		# excess in a following time step)

