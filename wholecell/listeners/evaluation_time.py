#!/usr/bin/env python

"""
EvaluationTime

Large-scale, low-overhead evaluation time tracker for process/state operations.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

class EvaluationTime(wholecell.listeners.listener.Listener):
	""" EvaluationTime """

	_name = 'EvaluationTime'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EvaluationTime, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(EvaluationTime, self).initialize(sim, kb)

		self.stateNames = sim.states.keys()
		self.processNames = sim.processes.keys()

		self.nStates = len(sim.states)
		self.nProcesses = len(sim.processes)

		# State evaluation times
		self.updateQueries_times = None
		self.partition_times = None
		self.merge_times = None

		self.updateQueries_total = None
		self.partition_total = None
		self.merge_total = None

		# Process evaluation times
		self.calculateRequest_times = None
		self.evolveState_times = None

		self.calculateRequest_total = None
		self.evolveState_total = None


	# Allocate memory
	def allocate(self):
		super(EvaluationTime, self).allocate()

		# State evaluation times
		self.updateQueries_times = np.zeros(self.nStates, np.float64)
		self.partition_times = np.zeros_like(self.updateQueries_times)
		self.merge_times = np.zeros_like(self.updateQueries_times)

		self.updateQueries_total = 0
		self.partition_total = 0
		self.merge_total = 0

		# Process evaluation times
		self.calculateRequest_times = np.zeros(self.nProcesses, np.float64)
		self.evolveState_times = np.zeros_like(self.calculateRequest_times)

		self.calculateRequest_total = 0
		self.evolveState_total = 0


	def update(self):
		self.updateQueries_total = self.updateQueries_times.sum()
		self.partition_total = self.partition_times.sum()
		self.merge_total = self.merge_times.sum()

		self.calculateRequest_total = self.calculateRequest_times.sum()
		self.evolveState_total = self.evolveState_times.sum()
		

	def pytablesCreate(self, h5file, expectedRows):

		# Columns
		dtype = {
			"time": tables.Int64Col(),
			"updateQueries_times": tables.Float64Col(self.nStates),
			"partition_times": tables.Float64Col(self.nStates),
			"merge_times": tables.Float64Col(self.nStates),
			"calculateRequest_times": tables.Float64Col(self.nProcesses),
			"evolveState_times": tables.Float64Col(self.nProcesses),
			"updateQueries_total": tables.Float64Col(),
			"partition_total": tables.Float64Col(),
			"merge_total": tables.Float64Col(),
			"calculateRequest_total": tables.Float64Col(),
			"evolveState_total": tables.Float64Col(),
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

		table.attrs.stateNames = self.stateNames
		table.attrs.processNames = self.processNames


	def pytablesAppend(self, h5file):

		table = h5file.get_node("/", self._name)
		entry = table.row

		entry["time"] = self.timeStep()
		entry["updateQueries_times"] = self.updateQueries_times
		entry["partition_times"] = self.partition_times
		entry["merge_times"] = self.merge_times
		entry["calculateRequest_times"] = self.calculateRequest_times
		entry["evolveState_times"] = self.evolveState_times
		entry["updateQueries_total"] = self.updateQueries_total
		entry["partition_total"] = self.partition_total
		entry["merge_total"] = self.merge_total
		entry["calculateRequest_total"] = self.calculateRequest_total
		entry["evolveState_total"] = self.evolveState_total

		entry.append()

		table.flush()
