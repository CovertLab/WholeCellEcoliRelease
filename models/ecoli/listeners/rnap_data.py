#!/usr/bin/env python

"""
RnapData

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/15
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class RnapData(wholecell.listeners.listener.Listener):
	""" RnapData """

	_name = 'RnapData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnapData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(RnapData, self).initialize(sim, kb)

		# Computed, saved attributes
		self.stallingRateTotal = None
		self.stallingRateMean = None
		self.stallingRateStd = None
		self.fractionStalled = None

		# Attributes broadcast by the processes
		self.rnapStallsFast = None
		self.rnapStallsSlow = None		
		self.ntpCountInSequence = None
		self.ntpCounts = None
		self.actualElongations = None
		self.expectedElongations = None
		self.nTerminatedFast = None
		self.nTerminatedSlow = None
		self.didTerminate = None
		self.didInitalize = None


		# Logged quantities
		self.registerLoggedQuantity(
			"Fraction\nrnaps\nstalled",
			"fractionStalled",
			".3f"
			)


	# Allocate memory
	def allocate(self):
		super(RnapData, self).allocate()

		# Computed, saved attributes
		self.stallingRateTotal = 0
		self.stallingRateMean = 0
		self.stallingRateStd = 0
		self.fractionStalled = 0

		# Attributes broadcast by the PolypeptideElongation process
		self.rnapStallsFast = np.zeros(0, np.int64)
		self.rnapStallsSlow = np.zeros(0, np.int64)
		self.ntpCountInSequence = np.zeros(21, np.int64)
		self.ntpCounts = np.zeros(21, np.int64)
		self.actualElongations = 0
		self.expectedElongations = 0
		self.nTerminatedFast = 0
		self.nTerminatedSlow = 0
		self.didTerminate = 0
		self.didInitalize = 0

	def update(self):
		if self.rnapStallsFast.size:
			# TODO: divide rates by time step length
			self.stallingRateFastTotal = self.rnapStallsFast.sum()
			self.stallingRateFastMean = self.rnapStallsFast.mean()
			self.stallingRateFastStd = self.rnapStallsFast.std()
			self.fractionStalledFast = (self.rnapStallsFast > 0).mean()

		else:
			self.stallingRateFastTotal = 0
			self.stallingRateFastMean = 0
			self.stallingRateFastStd = 0
			self.fractionStalledFast = 0

		if self.rnapStallsSlow.size:
			# TODO: divide rates by time step length
			self.stallingRateSlowTotal = self.rnapStallsSlow.sum()
			self.stallingRateSlowMean = self.rnapStallsSlow.mean()
			self.stallingRateSlowStd = self.rnapStallsSlow.std()
			self.fractionStalledSlow = (self.rnapStallsSlow > 0).mean()

		else:
			self.stallingRateSlowTotal = 0
			self.stallingRateSlowMean = 0
			self.stallingRateSlowStd = 0
			self.fractionStalledSlow = 0


	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			stallingRateFastTotal = self.stallingRateFastTotal,
			stallingRateFastMean = self.stallingRateFastMean,
			stallingRateFastStd = self.stallingRateFastStd,
			fractionStalledFast = self.fractionStalledFast,
			stallingRateSlowTotal = self.stallingRateSlowTotal,
			stallingRateSlowMean = self.stallingRateSlowMean,
			stallingRateSlowStd = self.stallingRateSlowStd,
			fractionStalledSlow = self.fractionStalledFast,
			ntpCountInSequence = self.ntpCountInSequence,
			ntpCounts = self.ntpCounts,
			actualElongations = self.actualElongations,
			expectedElongations = self.expectedElongations,
			nTerminatedFast = self.nTerminatedFast,
			nTerminatedSlow = self.nTerminatedSlow,
			didTerminate = self.didTerminate,
			didInitalize = self.didInitalize,
			)