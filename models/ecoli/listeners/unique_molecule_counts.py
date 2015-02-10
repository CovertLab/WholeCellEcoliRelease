#!/usr/bin/env python

"""
UniqueMoleculeCounts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

# TODO: move to the wholecell package & write interface such that it will
# function without requiring the state (will save an empty file)

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class UniqueMoleculeCounts(wholecell.listeners.listener.Listener):
	""" UniqueMoleculeCounts """

	_name = 'UniqueMoleculeCounts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(UniqueMoleculeCounts, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniqueMoleculeCounts, self).initialize(sim, kb)

		self.uniqueMolecules = sim.states['UniqueMolecules']

		self.uniqueMoleculeCounts = None


	# Allocate memory
	def allocate(self):
		super(UniqueMoleculeCounts, self).allocate()

		# TODO: add interface to unique objects container
		self.uniqueMoleculeCounts = np.zeros(
			len(self.uniqueMolecules.container._names),
			np.int64
			)


	def update(self):
		# TODO: add interface to unique objects container

		for i in xrange(self.uniqueMoleculeCounts.size):
			self.uniqueMoleculeCounts[i] = (
				self.uniqueMolecules.container._collections[i]["_entryState"]
				== self.uniqueMolecules.container._entryActive
				).sum()


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			uniqueMoleculeIds = self.uniqueMolecules.container._names
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			uniqueMoleculeCounts = self.uniqueMoleculeCounts,
			)
