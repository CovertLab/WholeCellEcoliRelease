#!/usr/bin/env python

"""
EquilibriumListener

Records dynamics of equilibrium output.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class EquilibriumListener(wholecell.listeners.listener.Listener):
	""" EquilibriumListener """

	_name = "EquilibriumListener"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EquilibriumListener, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(EquilibriumListener, self).initialize(sim, sim_data)

		self.complexIDs = sim_data.process.equilibrium.ids_complexes
		self.reactionIDs = sim_data.process.equilibrium.rxnIds


	# Allocate memory
	def allocate(self):
		super(EquilibriumListener, self).allocate()

		self.reactionRates = np.zeros(len(self.reactionIDs), np.float64)


	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionRates': 'reactionIDs'}

		tableWriter.writeAttributes(
			complexIDs = self.complexIDs,
			reactionIDs = self.reactionIDs,
			subcolumns = subcolumns)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionRates = self.reactionRates,
			)
