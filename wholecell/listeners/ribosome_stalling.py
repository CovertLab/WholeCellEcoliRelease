#!/usr/bin/env python

"""
RibosomeStalling

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

class RibosomeStalling(wholecell.listeners.listener.Listener):
	""" RibosomeStalling """

	_name = 'RibosomeStalling'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RibosomeStalling, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(RibosomeStalling, self).initialize(sim, kb)

		# Computed, saved attributes
		self.stallingRateTotal = None
		self.stallingRateMean = None
		self.stallingRateStd = None
		self.fractionStalled = None

		# Attributes broadcast by the UniquePolypeptideElongation process
		self.ribosomeStalls = None


	# Allocate memory
	def allocate(self):
		super(RibosomeStalling, self).allocate()

		# Computed, saved attributes
		self.stallingRateTotal = np.nan
		self.stallingRateMean = np.nan
		self.stallingRateStd = np.nan
		self.fractionStalled = np.nan

		# Attributes broadcast by the UniquePolypeptideElongation process
		self.ribosomeStalls = np.zeros(0, np.int64)


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
			"timeStep": tables.Int64Col(),
			"stallingRateTotal": tables.Float64Col(),
			"stallingRateMean": tables.Float64Col(),
			"stallingRateStd": tables.Float64Col(),
			"fractionStalled": tables.Float64Col(),
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

		entry["timeStep"] = self.timeStep()
		entry["stallingRateTotal"] = self.stallingRateTotal
		entry["stallingRateMean"] = self.stallingRateMean
		entry["stallingRateStd"] = self.stallingRateStd
		entry["fractionStalled"] = self.fractionStalled

		entry.append()

		table.flush()

	# TODO: load method
