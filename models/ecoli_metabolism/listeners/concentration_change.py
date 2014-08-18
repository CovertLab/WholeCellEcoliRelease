#!/usr/bin/env python

"""
ConcentrationChange

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/18/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

class ConcentrationChange(wholecell.listeners.listener.Listener):
	""" ConcentrationChange """

	_name = "ConcentrationChange"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ConcentrationChange, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(ConcentrationChange, self).initialize(sim, kb)

		self.metabolism = sim.processes["Metabolism"]

		self.moleculeIDs = kb.metabolitePoolIDs


	# Allocate memory
	def allocate(self):
		super(ConcentrationChange, self).allocate()
		
		self.concentrationChange = np.zeros(len(self.moleculeIDs), np.float64)


	def pytablesCreate(self, h5file, expectedRows):

		# Columns
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"concentrationChange": tables.Float64Col(self.concentrationChange.shape),
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
	
		groupNames = h5file.create_group(h5file.root,
			'names', 'Molecule names')

		h5file.create_array(groupNames, 'moleculeIDs', [str(s) for s in self.moleculeIDs])


	def pytablesAppend(self, h5file):

		table = h5file.get_node("/", self._name)
		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["concentrationChange"] = self.concentrationChange

		entry.append()

		table.flush()
