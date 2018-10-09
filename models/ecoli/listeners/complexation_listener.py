#!/usr/bin/env python

"""
ComplexationListener

Records dynamics of complexation output.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class ComplexationListener(wholecell.listeners.listener.Listener):
	""" ComplexationListener """

	_name = "ComplexationListener"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ComplexationListener, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ComplexationListener, self).initialize(sim, sim_data)

		self.complexIDs = sim_data.process.complexation.ids_complexes
		self.reactionIDs = sim_data.process.complexation.ids_reactions


	# Allocate memory
	def allocate(self):
		super(ComplexationListener, self).allocate()

		self.reactionRates = np.zeros(len(self.reactionIDs), np.float64)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			complexIDs = self.complexIDs,
			reactionIDs = self.reactionIDs,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionRates = self.reactionRates,
			)
