#!/usr/bin/env python

"""
FreeProduction

Provides the simulation with resources up to some desired level, as a function
of some initial count and the doubling time.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/5/2014
"""

import numpy

import wholecell.sim.process.Process

class FreeProduction(wholecell.sim.process.Process.Process):
	""" FreeProduction """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "FreeProduction",
			"name": "FreeProduction",
		}

		self.molIDs = None
		self.initCounts = None
		self.time = None

		self.defaultInitCount = 1e6
		self.doublingTime = 1. * 3600

		super(FreeProduction, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(FreeProduction, self).initialize(sim, kb)

		freeMolecules = sim.freeMolecules if sim.freeMolecules is not None else ()

		self.molIDs = []
		self.initCounts = self.defaultInitCount * numpy.ones(len(freeMolecules))

		for i, (molID, initCount) in enumerate(freeMolecules):
			self.molIDs.append(molID)

			if initCount is not None:
				self.initCounts[i] = initCount

		mc = sim.states["MoleculeCounts"]

		self.mcPartition = mc.addPartition(self, self.molIDs, self.calcReq)
		self.mcView = mc.countsBulkViewNew(self.molIDs)

		self.time = sim.states['Time']


	def calcReq(self, request):
		pass


	# Calculate temporal evolution
	def evolveState(self):
		expectedCounts = self.initCounts * numpy.exp(numpy.log(2) / self.doublingTime * self.time.value)

		self.mcPartition.countsBulkIs(
			numpy.fmax(
				0,
				expectedCounts - self.mcView.countsBulk()
				)
			)
