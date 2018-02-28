#!/usr/bin/env python

"""
GrowthLimits

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/25/15
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class GrowthLimits(wholecell.listeners.listener.Listener):
	""" GrowthLimits """

	_name = 'GrowthLimits'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(GrowthLimits, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(GrowthLimits, self).initialize(sim, sim_data)

		# Computed, saved attributes
		self.aaIds = sim_data.moleculeGroups.aaIDs
		self.ntpIds = sim_data.moleculeGroups.ntpIds

	# Allocate memory
	def allocate(self):
		super(GrowthLimits, self).allocate()

		# For translation
		self.gtpPoolSize = 0
		self.gtpRequestSize = 0
		self.gtpAllocated = 0
		self.gtpUsed = 0

		self.activeRibosomeAllocated = 0

		self.aaPoolSize = np.zeros(len(self.aaIds), np.float64)
		self.aaRequestSize = np.zeros(len(self.aaIds), np.float64)
		self.aaAllocated = np.zeros(len(self.aaIds), np.float64)
		self.aasUsed = np.zeros(len(self.aaIds), np.float64)

		# For transcription
		self.ntpPoolSize = np.zeros(len(self.ntpIds), np.float64)
		self.ntpRequestSize = np.zeros(len(self.ntpIds), np.float64)
		self.ntpAllocated = np.zeros(len(self.ntpIds), np.float64)
		self.ntpUsed = np.zeros(len(self.ntpIds), np.float64)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		pass

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			gtpPoolSize = self.gtpPoolSize,
			gtpRequestSize = self.gtpRequestSize,
			gtpAllocated = self.gtpAllocated,
			gtpUsed = self.gtpUsed,
			activeRibosomeAllocated = self.activeRibosomeAllocated,
			aaPoolSize = self.aaPoolSize,
			aaRequestSize = self.aaRequestSize,
			aaAllocated = self.aaAllocated,
			aasUsed = self.aasUsed,
			ntpPoolSize = self.ntpPoolSize,
			ntpRequestSize = self.ntpRequestSize,
			ntpAllocated = self.ntpAllocated,
			ntpUsed = self.ntpUsed,
			)