#!/usr/bin/env python

"""
EvaluationTime

Large-scale, low-overhead evaluation time tracker for process/state operations.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

from __future__ import division

import time

import numpy as np

import wholecell.listeners.listener

class EvaluationTime(wholecell.listeners.listener.Listener):
	""" EvaluationTime """

	_name = 'EvaluationTime'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EvaluationTime, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(EvaluationTime, self).initialize(sim, sim_data)
		# Clock time
		self.clock_time = None

		self.state_names = sim.internal_states.keys()
		self.process_names = sim.processes.keys()
		self.listener_names = sim.listeners.keys()
		self.logger_names = sim.loggers.keys()

		self.n_states = len(sim.internal_states)
		self.n_processes = len(sim.processes)
		self.n_listeners = len(sim.listeners)
		self.n_loggers = len(sim.loggers)

		# State evaluation times
		self.update_queries_times = None
		self.partition_times = None
		self.merge_times = None
		self.calculate_mass_times = None
		self.update_queries_total = None
		self.partition_total = None
		self.merge_total = None
		self.calculate_mass_total = None

		# Process evaluation times
		self.calculate_request_times = None
		self.evolve_state_times = None
		self.calculate_request_total = None
		self.evolve_state_total = None

		# Listener evaluation times
		self.update_times = None
		self.update_total = None

		# Logger evaluation times
		self.append_times = None
		self.append_total = None


	# Allocate memory
	def allocate(self):
		super(EvaluationTime, self).allocate()

		# Clock time
		self.clock_time = 0

		# State evaluation times
		self.update_queries_times = np.zeros(self.n_states, np.float64)
		self.partition_times = np.zeros_like(self.update_queries_times)
		self.merge_times = np.zeros_like(self.update_queries_times)
		self.calculate_mass_times = np.zeros_like(self.update_queries_times)
		self.update_queries_total = 0
		self.partition_total = 0
		self.merge_total = 0
		self.calculate_mass_total = 0

		# Process evaluation times
		self.calculate_request_times = np.zeros(self.n_processes, np.float64)
		self.evolve_state_times = np.zeros_like(self.calculate_request_times)
		self.calculate_request_total = 0
		self.evolve_state_total = 0

		# Listener evaluation times
		self.update_times = np.zeros(self.n_listeners, np.float64)
		self.update_total = 0

		# Logger evaluation times
		self.append_times = np.zeros(self.n_loggers, np.float64)
		self.append_total = 0


	def reset_evaluation_times(self):
		"""
		Reset evaluation time vectors of process and state methods that are
		run for each tier in the process hierarchy to zero.
		"""
		self.update_queries_times.fill(0)
		self.partition_times.fill(0)
		self.merge_times.fill(0)

		self.calculate_request_times.fill(0)
		self.evolve_state_times.fill(0)


	def update(self):
		self.clock_time = time.time()

		self.update_queries_total = self.update_queries_times.sum()
		self.partition_total = self.partition_times.sum()
		self.merge_total = self.merge_times.sum()
		self.calculate_mass_total = self.calculate_mass_times.sum()

		self.calculate_request_total = self.calculate_request_times.sum()
		self.evolve_state_total = self.evolve_state_times.sum()

		self.update_total = self.update_times.sum()

		self.append_total = self.append_times.sum()


	def tableCreate(self, tableWriter):
		# Handle the edge case of a simulation with no processes
		if self.n_processes == 0:
			return

		tableWriter.writeAttributes(
			state_names = self.state_names,
			process_names = self.process_names,
			listener_names = self.listener_names,
			logger_names = self.logger_names,
			)


	def tableAppend(self, tableWriter):
		# Handle the edge case of a simulation with no processes
		if self.n_processes == 0:
			return

		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			clock_time=self.clock_time,
			update_queries_times = self.update_queries_times,
			partition_times = self.partition_times,
			merge_times = self.merge_times,
			calculate_mass_times = self.calculate_mass_times,
			calculate_request_times = self.calculate_request_times,
			evolve_state_times = self.evolve_state_times,
			update_times = self.update_times,
			append_times = self.append_times,
			update_queries_total = self.update_queries_total,
			partition_total = self.partition_total,
			merge_total = self.merge_total,
			calculate_mass_total = self.calculate_mass_total,
			calculate_request_total = self.calculate_request_total,
			evolve_state_total = self.evolve_state_total,
			update_total = self.update_total,
			append_total = self.append_total,
			)
