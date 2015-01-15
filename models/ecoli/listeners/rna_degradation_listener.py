#!/usr/bin/env python

"""
RnaDegradationListener

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/15
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class RnaDegradationListener(wholecell.listeners.listener.Listener):
	""" RnaDegradationListener """

	_name = 'RnaDegradationListener'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnaDegradationListener, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradationListener, self).initialize(sim, kb)

		# Computed, saved attributes

		# Attributes broadcast by the RnaDegradation process
		self.countRnaDegraded = None
		self.nucleotidesFromDegradation = None

		# Logged quantities
		# self.registerLoggedQuantity(
		# 	"Fraction\nribosomes\nstalled",
		# 	"fractionStalled",
		# 	".3f"
		# 	)


	# Allocate memory
	def allocate(self):
		super(RnaDegradationListener, self).allocate()

		# Computed, saved attributes

		# Attributes broadcast by the RnaDegradation process
		self.countRnaDegraded = np.nan
		self.nucleotidesFromDegradation = np.nan

	def update(self):
		pass
		# if self.ribosomeStalls.size:
		# 	# TODO: divide rates by time step length
		# 	self.stallingRateTotal = self.ribosomeStalls.sum()
		# 	self.stallingRateMean = self.ribosomeStalls.mean()
		# 	self.stallingRateStd = self.ribosomeStalls.std()
		# 	self.fractionStalled = (self.ribosomeStalls > 0).mean()

		# else:
		# 	self.stallingRateTotal = np.nan
		# 	self.stallingRateMean = np.nan
		# 	self.stallingRateStd = np.nan
		# 	self.fractionStalled = np.nan


	def pytablesCreate(self, h5file, expectedRows):
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"countRnaDegraded": tables.Float64Col(),
			"nucleotidesFromDegradation": tables.Float64Col(),
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
		entry["countRnaDegraded"] = self.countRnaDegraded
		entry["nucleotidesFromDegradation"] = self.nucleotidesFromDegradation

		entry.append()

		table.flush()

	# TODO: load method
