#!/usr/bin/env python

"""
EnzymeKinetics

EnzymeKinetics listener. Tracks information about enzyme kinetics.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/02/2015
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

class EnzymeKinetics(wholecell.listeners.listener.Listener):
	""" EnzymeKinetics """

	_name = 'EnzymeKinetics'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EnzymeKinetics, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, kb):
		super(EnzymeKinetics, self).initialize(sim, kb)

		self.metabolism = sim.processes["Metabolism"]

	# Allocate memory
	# In case things are of unknown size, write them here
	# Dummy values for what will be writen to output table
	# prep variables with zeros or NaNs, with correct size
	# to be filled later
	def allocate(self):
		super(EnzymeKinetics, self).allocate()

		self.reactionRates = np.zeros(len(self.metabolism.fba.reactionIDs()), np.float64)
		self.reactionIDs = self.metabolism.fba.reactionIDs()

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			# metaboliteIds = self.metaboliteIds,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			reactionIDs = self.reactionIDs,
			reactionRates = self.reactionRates,
			)