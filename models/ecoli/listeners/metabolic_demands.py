#!/usr/bin/env python

"""
MetabolicDemands

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/10/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

class MetabolicDemands(wholecell.listeners.listener.Listener):
	""" MetabolicDemands """

	_name = "MetabolicDemands"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(MetabolicDemands, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(MetabolicDemands, self).initialize(sim, kb)

		self.bulkMolecules = sim.states["BulkMolecules"]

		self.processes = sim.processes

		self.cellCycleLen = kb.cellCycleLen.to('s').magnitude

		self.metaboliteIds = kb.wildtypeBiomass["metaboliteId"]

		moleculeIds = kb.bulkMolecules["moleculeId"]

		self.metaboliteIndexes = np.array([
			np.where(moleculeIds == metaboliteId)[0][0]
			for metaboliteId in self.metaboliteIds
			])

		self.metaboliteRequests = None


	# Allocate memory
	def allocate(self):
		super(MetabolicDemands, self).allocate()

		self.metaboliteRequests = np.zeros(
			(self.metaboliteIndexes.size, len(self.processes)),
			np.int64
			)


	def updatePostRequest(self):
		self.metaboliteRequests = self.bulkMolecules._countsRequested[
			self.metaboliteIndexes, :]


	def pytablesCreate(self, h5file, expectedRows):

		shape = self.metaboliteRequests.shape

		# Columns
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"metaboliteRequests": tables.UInt64Col(shape),
			}

		# Create table
		table = h5file.create_table(
			h5file.root,
			self._name,
			dtype,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		# Store units as metadata
		table.attrs.metaboliteIds = self.metaboliteIds
		table.attrs.processIds = self.processes.keys()


	def pytablesAppend(self, h5file):
		table = h5file.get_node("/", self._name)
		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["metaboliteRequests"] = self.metaboliteRequests

		entry.append()

		table.flush()
