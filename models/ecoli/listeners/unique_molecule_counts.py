"""
UniqueMoleculeCounts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

# TODO: move to the wholecell package & write interface such that it will
# function without requiring the state (will save an empty file)

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class UniqueMoleculeCounts(wholecell.listeners.listener.Listener):
	""" UniqueMoleculeCounts """

	_name = 'UniqueMoleculeCounts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(UniqueMoleculeCounts, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(UniqueMoleculeCounts, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']


	# Allocate memory
	def allocate(self):
		super(UniqueMoleculeCounts, self).allocate()

		self.uniqueMoleculeCounts = np.zeros(
			len(self.uniqueMolecules.container.objectNames()),
			np.int64
			)


	def update(self):
		self.uniqueMoleculeCounts = self.uniqueMolecules.container.counts()


	def tableCreate(self, tableWriter):
		objectNames = self.uniqueMolecules.container.objectNames()
		subcolumns = {
			'uniqueMoleculeCounts': 'objectNames'}

		tableWriter.writeAttributes(
			uniqueMoleculeIds = self.uniqueMolecules.container.objectNames(),
			objectNames = objectNames,
			subcolumns = subcolumns)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			uniqueMoleculeCounts = self.uniqueMoleculeCounts,
			)
