#!/usr/bin/env python

"""
UniqueMoleculeCounts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

class UniqueMoleculeCounts(wholecell.listeners.listener.Listener):
	""" UniqueMoleculeCounts """

	_name = 'UniqueMoleculeCounts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(UniqueMoleculeCounts, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniqueMoleculeCounts, self).initialize(sim, kb)

		self.uniqueMolecules = sim.states['UniqueMolecules']

		self.uniqueMoleculeCounts = None



	# Allocate memory
	def allocate(self):
		super(UniqueMoleculeCounts, self).allocate()

		# TODO: add interface to unique objects container
		self.uniqueMoleculeCounts = np.zeros(
			len(self.uniqueMolecules.container._names),
			np.int64
			)


	def update(self):
		# TODO: add interface to unique objects container

		for i in xrange(self.uniqueMoleculeCounts.size):
			self.uniqueMoleculeCounts[i] = (
				self.uniqueMolecules.container._collections[i]["_entryState"]
				== self.uniqueMolecules.container._entryActive
				).sum()



	def pytablesCreate(self, h5file, expectedRows):
		size = self.uniqueMoleculeCounts.size
		# Columns
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"uniqueMoleculeCounts": tables.UInt64Col(size)
			}

		# Create table
		# TODO: Add compression options (using filters)
		table = h5file.create_table(
			h5file.root,
			self._name,
			dtype,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		table.attrs.uniqueMoleculeIds = self.uniqueMolecules.container._names


	def pytablesAppend(self, h5file):
		table = h5file.get_node("/", self._name)
		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["uniqueMoleculeCounts"] = self.uniqueMoleculeCounts

		entry.append()

		table.flush()
