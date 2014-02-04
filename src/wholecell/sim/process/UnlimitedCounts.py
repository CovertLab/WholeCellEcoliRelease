#!/usr/bin/env python

"""
UnlimitedCounts

A simplified submodel that maintains counts of molecules at some minimum level.
Used to test other submodels in isolation.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/4/2014
"""

import numpy

import wholecell.sim.process.Process

class UnlimitedCounts(wholecell.sim.process.Process.Process):
	""" UnlimitedCounts """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "UnlimitedCounts",
			"name": "UnlimitedCounts",
			"options": [
				"fixedMolecules",	# iterable of tuples of (ID, counts)
				]
			}

		# Options
		self.fixedMolecules = []

		# Partitions
		self.mcPartition = None

		# Counts
		self.molIDs = []
		self.minCounts = None

		# Constants
		self.defaultCount = 1e6

		super(UnlimitedCounts, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UnlimitedCounts, self).initialize(sim, kb)

		self.minCounts = self.defaultCount * numpy.ones(len(self.fixedMolecules), float)

		for i, (molID, count) in enumerate(self.fixedMolecules):
			self.molIDs.append(molID)

			if count is not None:
				self.minCounts[i] = count

		self.mcPartition = mc.addPartition(self, self.molIDs, self.calcRequest)
		self.mcView = mc.countsBulkViewNew(self.molIDs)


	def calcRequest(self, request):
		pass


	# Calculate temporal evolution
	def evolveState(self):
		self.mcPartition.countsBulkIs(
			numpy.fmax(0, self.minCounts - self.mcView.countsBulk())
			)
