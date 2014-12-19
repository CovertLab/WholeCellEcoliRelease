#!/usr/bin/env python

"""
RibosomeData

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/14
"""

from __future__ import division

import numpy as np
import tables

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

		# Attributes broadcast by the PolypeptideElongation process
		self.ribosomeStalls = None
		self.aaCountInSequence = None
		self.aaCounts = None
		self.trnasCapacity = None
		self.synthetaseCapacity = None
		self.actualElongations = None
		self.expectedElongations = None


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
		self.stallingRateTotal = np.nan
		self.stallingRateMean = np.nan
		self.stallingRateStd = np.nan
		self.fractionStalled = np.nan

		# Attributes broadcast by the PolypeptideElongation process
		self.ribosomeStalls = np.zeros(0, np.int64)
		self.aaCountInSequence = np.zeros(21, np.int64)
		self.aaCounts = np.zeros(21, np.int64)
		self.trnasCapacity = np.zeros(21, np.int64)
		self.synthetaseCapacity = np.zeros(21, np.int64)
		self.actualElongations = np.nan
		self.expectedElongations = np.nan

	def update(self):
		if self.ribosomeStalls.size:
			# TODO: divide rates by time step length
			self.stallingRateTotal = self.ribosomeStalls.sum()
			self.stallingRateMean = self.ribosomeStalls.mean()
			self.stallingRateStd = self.ribosomeStalls.std()
			self.fractionStalled = (self.ribosomeStalls > 0).mean()

		else:
			self.stallingRateTotal = np.nan
			self.stallingRateMean = np.nan
			self.stallingRateStd = np.nan
			self.fractionStalled = np.nan


	def pytablesCreate(self, h5file, expectedRows):
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"stallingRateTotal": tables.Float64Col(),
			"stallingRateMean": tables.Float64Col(),
			"stallingRateStd": tables.Float64Col(),
			"fractionStalled": tables.Float64Col(),
			"aaCountInSequence": tables.Float64Col(self.aaCountInSequence.size),
			"aaCounts": tables.Float64Col(self.aaCounts.size),
			"trnasCapacity": tables.Float64Col(self.trnasCapacity.size),
			"synthetaseCapacity": tables.Float64Col(self.synthetaseCapacity.size),
			"actualElongations": tables.Float64Col(),
			"expectedElongations": tables.Float64Col(),
			}

		table = h5file.create_table(
			h5file.root,
			self._name,
			dtype,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			)


	def pytablesAppend(self, h5file):
		table = h5file.get_node("/", self._name)

		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["stallingRateTotal"] = self.stallingRateTotal
		entry["stallingRateMean"] = self.stallingRateMean
		entry["stallingRateStd"] = self.stallingRateStd
		entry["fractionStalled"] = self.fractionStalled
		entry["aaCountInSequence"] = self.aaCountInSequence
		entry["aaCounts"] = self.aaCounts
		entry["trnasCapacity"] = self.trnasCapacity
		entry["synthetaseCapacity"] = self.synthetaseCapacity
		entry["actualElongations"] = self.actualElongations
		entry["expectedElongations"] = self.expectedElongations

		entry.append()

		table.flush()

	# TODO: load method
