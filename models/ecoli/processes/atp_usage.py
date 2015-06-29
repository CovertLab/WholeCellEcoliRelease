#!/usr/bin/env python

"""
AtpUsage

Hydrolyze ATP (for non-growth associated maintenance)

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
from wholecell.utils import units

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
		self.nAvogadro = kb.constants.nAvogadro.asNumber(1 / units.mol)

		moleculeIds = ["ATP[c]", "WATER[c]", "Pi[c]", "ADP[c]", "PROTON[c]"]
		self.molecules = self.bulkMoleculesView(moleculeIds)
		self.reactants = self.bulkMoleculesView(["ATP[c]", "WATER[c]"])
		self.products = self.bulkMoleculesView(["ADP[c]", "Pi[c]", "PROTON[c]"])
		self.atp = self.bulkMoleculeView("ATP[c]")
		self.h2o = self.bulkMoleculeView("WATER[c]")
		# self.pi = self.bulkMoleculeView("Pi[c]")
		# self.adp = self.bulkMoleculeView("ADP[c]")
		# self.h = self.bulkMoleculeView("PROTON[c]")

		self.nongrowthAssociated_reactionsPerTimestep = (
			kb.constants.nonGrowthAssociatedMaintenance * kb.constants.nAvogadro
			).asNumber(1/units.fg/units.s) * self.timeStepSec

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_ATP_USAGE)


	def calculateRequest(self):
		mass = self.readFromListener("Mass", "dryMass")
		deltaMass = self.readFromListener("Mass", "growth")

		expectedReactions_nongrowthAssociated = (mass
			* self.nongrowthAssociated_reactionsPerTimestep)

		expectedReactions = expectedReactions_nongrowthAssociated

		self.reactants.requestIs(
			stochasticRound(self.randomState, expectedReactions)
			)


	def evolveState(self):
		mass = self.readFromListener("Mass", "dryMass")
		deltaMass = self.readFromListener("Mass", "growth")

		expectedReactions_nongrowthAssociated = (mass
			* self.nongrowthAssociated_reactionsPerTimestep)

		expectedReactions = expectedReactions_nongrowthAssociated

		atpsHydrolyzed = np.fmin(
			self.atp.count(),
			self.h2o.count()
			)

		self.reactants.countsDec(atpsHydrolyzed)
		self.products.countsInc(atpsHydrolyzed)

		if atpsHydrolyzed < np.floor(expectedReactions):
			warnings.warn("Maintenance not satisfied; the cell may be growing too fast")

			print "did not meet maintenance ({})".format(expectedReactions - atpsHydrolyzed)

			# TODO: flag some sort of listener that tracks phenomenological
			# observations like cell death/division

		# TODO: consider implementing (N)GAM as a metabolism constraint
		# TODO: let (N)GAM roll over to the next time step, just in case
		# there is some stochasticity (i.e. ATP limiting in one time step, in
		# excess in a following time step)

