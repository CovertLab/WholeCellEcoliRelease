from __future__ import absolute_import, division, print_function

import time
import random

from .inner import CellSimulation
from .outer import EnvironmentSimulation

class SimulationStub(CellSimulation):

	"""
	Provide a stub for the simulation.

	This class is intended to be a demonstration implementation of the CellSimulation
    interface from `agent/inner.py`.
	"""

	def __init__(self):
		self.local_time = 0
		self.local_set = False
		self.concentrations = {}
		self.environment_change = {}

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.concentrations = update['concentrations']

		# TODO(Ryan): consider if the outer agent should send a previous message
		#     establishing what keys are in this dict to avoid this preset state. 
		if not self.local_set:
			self.environment_change = {}
			for molecule in self.concentrations.iterkeys():
				self.environment_change[molecule] = 0
			self.local_set = True

	def run_incremental(self, run_until):
		# The following emulates simulations taking less time than the provided
		# `run_until` in order to validate that the environment as a whole proceeds
		# from time step to time step correctly.

		# (1 / interrupt_frequency) is how often the interrupt occurs
		interrupt_frequency = 6
		interrupt = random.randint(1, interrupt_frequency)
		interrupted = interrupt > interrupt_frequency - 1
		step = (run_until - self.local_time) / interrupt_frequency
		until = run_until
		if interrupted:
			until = random.randint(1, interrupt_frequency - 1) * step + self.local_time
		span = until - self.local_time

		# print('=== simulation | run_until: {}, until: {}, step: {}, span: {}'.format(run_until, until, step, span))

		time.sleep(random.randint(1, interrupt_frequency))

		for molecule in self.concentrations.iterkeys():
			self.environment_change[molecule] += random.randint(1, 6) * span

		self.local_time = until

	def generate_inner_update(self):
		return {
			'changes': self.environment_change}

	def finalize(self):
		pass

class EnvironmentStub(EnvironmentSimulation):

	"""
	Provide a stub for the environmental context.

	Like the stub above for simulations, this class is meant to demonstrate the
	EnvironmentSimulation interface defined in `agent/outer.py` that any environmental
	simulation must implement to fulfill the conditions of being a coordinating
	external context for integrating the changes of each simulation.

	Full interface described in `agent/outer.py`.
	"""

	def __init__(self, volume, concentrations):
		self.global_time = 0
		self.volume = volume
		self.simulations = {}
		self.concentrations = concentrations
		self.run_for = 1.0

	def time(self):
		return self.global_time

	def add_simulation(self, agent_id, state):
		self.simulations[agent_id] = state

	def remove_simulation(self, agent_id):
		return self.simulations.pop(agent_id, {})

	def apply_inner_update(self, update, now):
		self.simulations.update(update)

		for agent_id, simulation in self.simulations.iteritems():
			if simulation['time'] <= now:
				state = simulation['state']
				for molecule, change in state['changes'].iteritems():
					self.concentrations[molecule] += change

	def run_for_time(self):
		return self.run_for

	def run_incremental(self, run_until):
		time.sleep(2)
		self.global_time = run_until

	def generate_outer_update(self, now):
		state = {}
		for agent_id, simulation in self.simulations.iteritems():
			if simulation['time'] <= now:
				state[agent_id] = {}
				state[agent_id]['concentrations'] = self.concentrations

		return state

