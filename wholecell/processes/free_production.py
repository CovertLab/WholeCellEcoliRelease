#!/usr/bin/env python

"""
FreeProduction

Provides the simulation with resources up to some desired level, as a function
of some initial count and the doubling time.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/5/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class FreeProduction(wholecell.processes.process.Process):
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
		self.doublingTime = 1. * 3600 # TOKB

		super(FreeProduction, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(FreeProduction, self).initialize(sim, kb)

		freeMolecules = sim.freeMolecules if sim.freeMolecules is not None else ()

		self.molIDs = []
		self.initCounts = self.defaultInitCount * np.ones(len(freeMolecules))

		for i, (molID, initCount) in enumerate(freeMolecules):
			self.molIDs.append(molID)

			if initCount is not None:
				self.initCounts[i] = initCount

		mc = sim.states["BulkMolecules"]

		self.bulkMoleculesPartition.initialize(self.molIDs)
		self.mcView = mc.countsView(self.molIDs)

		self.time = sim.states['Time']

		# Views
		self.molecules = self.bulkMoleculesView(self.molIDs)


	def requestMoleculeCounts(self):
		# No request, since it only produces molecules
		pass


	def calculateRequest(self):
		# No request, since it only produces molecules
		pass


	# Calculate temporal evolution
	def evolveState(self):
		expectedCounts = self.initCounts * np.exp(np.log(2) / self.doublingTime * self.time.value)

		self.bulkMoleculesPartition.countsIs(
			np.fmax(
				0,
				expectedCounts - self.mcView.counts()
				)
			)
