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
		self.rnapStalls = None		
		self.ntpCountInSequence = None
		self.ntpCounts = None
		self.actualElongations = None
		self.expectedElongations = None
		self.nTerminated = None
		self.didTerminate = None
		self.didInitialize = None
		self.terminationLoss = None


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
		self.rnapStalls = np.zeros(0, np.int64)
		self.ntpCountInSequence = np.zeros(21, np.int64)
		self.ntpCounts = np.zeros(21, np.int64)
		self.actualElongations = 0
		self.expectedElongations = 0
		self.nTerminated = 0
		self.didTerminate = 0
		self.didInitialize = 0
		self.terminationLoss = 0

	def update(self):
		if self.rnapStalls.size:
			# TODO: divide rates by time step length
			self.stallingRateTotal = self.rnapStalls.sum()
			self.stallingRateMean = self.rnapStalls.mean()
			self.stallingRateStd = self.rnapStalls.std()
			self.fractionStalled = (self.rnapStalls > 0).mean()

		else:
			self.stallingRateTotal = 0
			self.stallingRateMean = 0
			self.stallingRateStd = 0
			self.fractionStalled = 0


	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			stallingRateTotal = self.stallingRateTotal,
			stallingRateMean = self.stallingRateMean,
			stallingRateStd = self.stallingRateStd,
			fractionStalled = self.fractionStalled,
			ntpCountInSequence = self.ntpCountInSequence,
			ntpCounts = self.ntpCounts,
			actualElongations = self.actualElongations,
			expectedElongations = self.expectedElongations,
			nTerminated = self.nTerminated,
			didTerminate = self.didTerminate,
			didInitialize = self.didInitialize,
			terminationLoss = self.terminationLoss
			)