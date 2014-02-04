#!/usr/bin/env python

"""
GrowingCounts

A simplified submodel that maintains the counts of molecules at some level, 
growing exponentially in time with the expected cellular growth rate.  Used to
test other submodels in isolation.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/4/2014
"""

import numpy

import wholecell.sim.process.Process

class GrowingCounts(wholecell.sim.process.Process.Process):
	""" GrowingCounts """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "GrowingCounts",
			"name": "GrowingCounts",
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
		self.initialCounts = None

		# Constants
		self.defaultCount = 1e6

		super(UnlimitedCounts, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UnlimitedCounts, self).initialize(sim, kb)

		self.time = sim.states['Time']

		self.initialCounts = self.defaultCount * numpy.ones(len(self.fixedMolecules), float)

		for i, (molID, count) in enumerate(self.fixedMolecules):
			self.molIDs.append(molID)
			self.initialCounts[i] = count

		self.mcPartition = mc.addPartition(self, self.molIDs, self.calcRequest)
		self.mcView = mc.countsBulkViewNew(self.molIDs)


	def calcRequest(self, request):
		pass


	# Calculate temporal evolution
	def evolveState(self):
		expectedCounts = numpy.exp(
			numpy.log(2)/self.cellCycleLen*self.time.value
			) * self.initialCounts

		self.mcPartition.countsBulkIs(
			numpy.fmax(0, expectedCounts - self.mcView.countsBulk())
			)
