#!/usr/bin/env python

"""
RibosomeData

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/14
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class RibosomeData(wholecell.listeners.listener.Listener):
	""" RibosomeData """

	_name = 'RibosomeData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RibosomeData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(RibosomeData, self).initialize(sim, kb)

		# Computed, saved attributes
		self.stallingRateTotal = None
		self.stallingRateMean = None
		self.stallingRateStd = None
		self.fractionStalled = None

		# Attributes broadcast by the processes
		self.ribosomeStalls = None
		self.aaCountInSequence = None
		self.aaCounts = None
		self.actualElongations = None
		self.expectedElongations = None
		self.didTerminate = None
		self.didInitalize = None
		self.terminationLoss = None

		# Logged quantities
		self.registerLoggedQuantity(
			"Fraction\nribosomes\nstalled",
			"fractionStalled",
			".3f"
			)


	# Allocate memory
	def allocate(self):
		super(RibosomeData, self).allocate()

		# Computed, saved attributes
		self.stallingRateTotal = 0
		self.stallingRateMean = 0
		self.stallingRateStd = 0
		self.fractionStalled = 0

		# Attributes broadcast by the PolypeptideElongation process
		self.ribosomeStalls = np.zeros(0, np.int64)
		self.aaCountInSequence = np.zeros(21, np.int64)
		self.aaCounts = np.zeros(21, np.int64)
		self.actualElongations = 0
		self.expectedElongations = 0
		self.didTerminate = 0
		self.didInitalize = 0
		self.terminationLoss = 0

	def update(self):
		if self.ribosomeStalls.size:
			# TODO: divide rates by time step length
			self.stallingRateTotal = self.ribosomeStalls.sum()
			self.stallingRateMean = self.ribosomeStalls.mean()
			self.stallingRateStd = self.ribosomeStalls.std()
			self.fractionStalled = (self.ribosomeStalls > 0).mean()

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
			aaCountInSequence = self.aaCountInSequence,
			aaCounts = self.aaCounts,
			actualElongations = self.actualElongations,
			expectedElongations = self.expectedElongations,
			didTerminate = self.didTerminate,
			didInitalize = self.didInitalize,
			terminationLoss = self.terminationLoss,
			)
