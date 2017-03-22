#!/usr/bin/env python

"""
AtpUsage

Hydrolyze ATP (for non-growth associated maintenance)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014

TODO: consider implementing (N)GAM as a metabolism constraint
TODO: let (N)GAM roll over to the next time step, just in case there is some stochasticity
(i.e. ATP limiting in one time step, in excess in a following time step)
TODO: Flag listner that tracks cell death when maintenance is not met
"""

from __future__ import division

import warnings

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_ATP_USAGE
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

VERBOSE = False

class AtpUsage(wholecell.processes.process.Process):
	""" AtpUsage """

	_name = "AtpUsage"

	def __init__(self):
		super(AtpUsage, self).__init__()

	def initialize(self, sim, sim_data):
		super(AtpUsage, self).initialize(sim, sim_data)

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)

		# Create Views on necessary metabolites for ATP hydrolysis reaction
		# ATP + H2O --> ADP + Pi + H
		moleculeIds = ["ATP[c]", "WATER[c]", "PI[c]", "ADP[c]", "PROTON[c]"]
		self.molecules = self.bulkMoleculesView(moleculeIds)
		self.reactants = self.bulkMoleculesView(["ATP[c]", "WATER[c]"])
		self.products = self.bulkMoleculesView(["ADP[c]", "PI[c]", "PROTON[c]"])
		self.atp = self.bulkMoleculeView("ATP[c]")
		self.h2o = self.bulkMoleculeView("WATER[c]")

		# Load parameter for non-growth associated maintenance from SimulationData
		self.nonGrowthMaintenance = (sim_data.constants.nonGrowthAssociatedMaintenance * sim_data.constants.nAvogadro
			).asNumber(1/units.fg/units.s)

		# Set request priority for state partitioning. This ensures that ATP hydrolysis for
		# cell maintenance always has a priority request and maintenance requirements are met
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_ATP_USAGE)


	def calculateRequest(self):
		# Get dry mass of the cell and growth of the cell in last simulation step
		mass = self.readFromListener("Mass", "dryMass")
		deltaMass = self.readFromListener("Mass", "growth")

		# NGAM is in units of per cell mass and per time this converts this to reaction counts by multiplying
		# by cell mass and the time step for the current simulation step. This can be a fraction.
		expectedReactions = (mass
			* self._nongrowthAssociated_reactionsPerTimestep())

		# Rounds expectedReactions stochastically to the nearest integer value (weighted by decimal)
		# and requests correct number of reactant metabolites to carry out this number of reactions
		self.reactants.requestIs(
			stochasticRound(self.randomState, expectedReactions)
			)

	def evolveState(self):
		# Get dry mass of the cell and growth of the cell in last simulation step
		mass = self.readFromListener("Mass", "dryMass")
		deltaMass = self.readFromListener("Mass", "growth")

		# NGAM is in units of per cell mass and per time this converts this to reaction counts by multiplying
		# by cell mass and the time step for the current simulation step. This can be a fraction.
		expectedReactions = (mass
			* self._nongrowthAssociated_reactionsPerTimestep())

		# Calculates how much ATP to hydrolyze limited by water counts
		atpsHydrolyzed = np.fmin(
			self.atp.count(),
			self.h2o.count()
			)

		# Logs data
		if VERBOSE: print "ATP hydrolyzed = %f" % atpsHydrolyzed
		self.writeToListener("ATPhydrolyzedUsageListener", "atpsHydrolyzed", atpsHydrolyzed)

		# Decrements reactant counts and increments product counts
		self.reactants.countsDec(atpsHydrolyzed)
		self.products.countsInc(atpsHydrolyzed)

		# If for some reason maintence is not met throw a warning
		if atpsHydrolyzed < np.floor(expectedReactions):
			warnings.warn("Maintenance not satisfied; the cell may be growing too fast")
			print "did not meet maintenance ({})".format(expectedReactions - atpsHydrolyzed)

	def _nongrowthAssociated_reactionsPerTimestep(self):
		return self.nonGrowthMaintenance * self.timeStepSec()