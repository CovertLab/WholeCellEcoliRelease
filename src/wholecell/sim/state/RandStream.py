#!/usr/bin/env python

"""
RandStream

Pseudorandom number generator state.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/31/2014
"""

import numpy
import tables

import wholecell.sim.state.State

class RandStream(wholecell.sim.state.State.State):
	""" RandStream """

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "RandStream",
			"name": "RandStream",
			"dynamics": ["value"]
			}

		self.value = None

		self.time = None

		super(RandStream, self).__init__(*args, **kwargs)

	def initialize(self, sim, kb):
		super(RandStream, self).initialize(sim, kb)

		self.time = sim.states["Time"]


	def calculate(self):
		self.value = self.randStream.state


	def pytablesCreate(self, h5file):
		stringLen = len(self.value[0])
		arrayShape = self.value[1].shape

		# Columns
		d = {
			"time": tables.Int64Col(),
			'string':tables.StringCol(stringLen), # string columns require an item length
			'integerArray':tables.Int32Col(arrayShape),
			'position':tables.UInt64Col(),
			'hasGauss':tables.UInt64Col(),
			'cachedGaussian':tables.Float64Col()
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self.meta["id"],
			d,
			title = self.meta["name"],
			filters = tables.Filters(complevel = 9, complib="zlib")
			)


	def pytablesAppend(self, h5file):
		[string, integerArray, position, hasGauss, cachedGaussian] = self.value # see numpy documentation for an explanation of these terms

		simTime = self.time.value

		t = h5file.get_node("/", self.meta["id"])
		entry = t.row

		entry["time"] = simTime
		entry['string'] = string
		entry['integerArray'] = integerArray
		entry['position'] = position
		entry['hasGauss'] = hasGauss
		entry['cachedGaussian'] = cachedGaussian

		entry.append()

		t.flush()


	def pytablesLoad(self, h5file, timePoint):
		entry = h5file.get_node('/', self.meta['id'])[timePoint]

		state = (entry['string'].read(), entry['integerArray'].read(),
			entry['position'].read(), entry['hasGauss'].read(),
			entry['cachedGaussian'].read())

		self.randStream.state = state
