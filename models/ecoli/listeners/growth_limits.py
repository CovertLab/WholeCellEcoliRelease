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

		# Attributes broadcast by the processes
		self.fractionAAsUsed = None
		self.fractionGtpLimit = None
		self.gtpPoolSize = None
		self.gtpRequestSize = None
		self.gtpAllocated = None
		self.gtpPerElongation = None


	# Allocate memory
	def allocate(self):
		super(GrowthLimits, self).allocate()

		self.fractionAAsUsed = 0.
		self.fractionGtpLimit = 0.
		self.gtpPoolSize = 0
		self.gtpRequestSize = 0
		self.gtpAllocated = 0
		self.gtpPerElongation = 0

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
			gtpPerElongation = self.gtpPerElongation,
			)