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
	def initialize(self, sim, kb):
		super(GrowthLimits, self).initialize(sim, kb)

		# Computed, saved attributes
		self.aaIds = kb.moleculeGroups.aaIDs

		# For transcription
		self.gtpPoolSize = None
		self.gtpRequestSize = None
		self.gtpAllocated = None
		self.gtpUsed = None

		self.aaPoolSize = None
		self.aaRequestSize = None
		self.aaAllocated = None
		self.aasUsed = None

	# Allocate memory
	def allocate(self):
		super(GrowthLimits, self).allocate()

		self.fractionAAsUsed = 0.
		self.fractionGtpLimit = 0.
		self.gtpPoolSize = 0
		self.gtpRequestSize = 0
		self.gtpAllocated = 0
		self.gtpUsed = 0

		self.aaPoolSize = np.zeros(len(self.aaIds), np.float64)
		self.aaRequestSize = np.zeros(len(self.aaIds), np.float64)
		self.aaAllocated = np.zeros(len(self.aaIds), np.float64)
		self.aasUsed = np.zeros(len(self.aaIds), np.float64)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		pass

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			fractionAAsUsed = self.fractionAAsUsed,
			fractionGtpLimit = self.fractionGtpLimit,
			gtpPoolSize = self.gtpPoolSize,
			gtpRequestSize = self.gtpRequestSize,
			gtpAllocated = self.gtpAllocated,
			gtpUsed = self.gtpUsed,
			aaPoolSize = self.aaPoolSize,
			aaRequestSize = self.aaRequestSize,
			aaAllocated = self.aaAllocated,
			aasUsed = self.aasUsed,
			)