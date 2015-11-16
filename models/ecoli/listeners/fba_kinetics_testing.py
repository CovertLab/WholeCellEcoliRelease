#!/usr/bin/env python

"""
FbaKineticsTesting

FbaKineticsTesting listener. Saves information to allow testing FBA with enzyme kinetics outside the main model.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/17/2015
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

class FbaKineticsTesting(wholecell.listeners.listener.Listener):
	""" FbaKineticsTesting """

	_name = 'FbaKineticsTesting'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(FbaKineticsTesting, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, kb):
		super(FbaKineticsTesting, self).initialize(sim, kb)

		self.metabolism = sim.processes["Metabolism"]

	# Allocate memory
	# In case things are of unknown size, write them here
	# Dummy values for what will be writen to output table
	# prep variables with zeros or NaNs, with correct size
	# to be filled later
	def allocate(self):
		super(FbaKineticsTesting, self).allocate()

		self.reactionRates = np.zeros(len(self.metabolism.fba.reactionIDs()), np.float64)
		self.perEnzymeRates = np.zeros(len(self.metabolism.fba.reactionIDs()), np.float64)
		self.reactionIDs = self.metabolism.fba.reactionIDs()

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = self.reactionIDs,

			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			reactionRates = self.reactionRates,
			perEnzymeRates = self.perEnzymeRates,
			)