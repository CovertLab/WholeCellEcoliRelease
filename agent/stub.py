from __future__ import absolute_import, division, print_function

import time
import numpy as np
import random

from .inner import CellSimulation
from .outer import EnvironmentSimulation

class SimulationStub(CellSimulation):

	"""
	Provide a stub for the simulation.

	This class is not meant to be functional but rather to implement in the simplest way 
	the interface the agents expect when interacting with the simulation.

	Full interface documented in `agent/inner.py`.
	"""

	def __init__(self):
		self.local_time = 0
		self.local_set = False
		self.concentrations = {}
		self.environment_change = {}

	def time(self):
		return self.local_time

	def initialize_local_environment(self):
		pass

	def set_local_environment(self, concentrations):
		self.concentrations = concentrations

		# TODO(Ryan): consider if the outer agent should send a previous message
		#     establishing what keys are in this dict to avoid this preset state. 
		if not self.local_set:
			self.environment_change = {}
			for molecule in self.concentrations.iterkeys():
				self.environment_change[molecule] = 0
			self.local_set = True

	def run_incremental(self, run_until):
		interrupt_frequency = 6
		interrupt = random.randint(1, interrupt_frequency)
		interrupted = interrupt > 5
		step = (run_until - self.local_time) / interrupt_frequency
		until = random.randint(
			1, interrupt_frequency - 1) * step + self.local_time if interrupted else run_until
		span = until - self.local_time

		print('================ simulation | run_until: {}, until: {}, step: {}, span: {}'.format(run_until, until, step, span))

		time.sleep(random.randint(1, interrupt_frequency))

		for molecule in self.concentrations.iterkeys():
			self.environment_change[molecule] += random.randint(1, 6) * span

		self.local_time = until

	def get_environment_change(self):
		return {
			'changes': self.environment_change,
			'time': self.time()}

	def finalize(self):
		pass

class EnvironmentStub(EnvironmentSimulation):

	"""
	Provide a stub for the environmental context.

	Like the stub above for simulations, this class is meant to demonstrate the interface
	any environmental simulation must have to fulfill the conditions of being a coordinating
	external context for integrating the changes of each simulation.

	Full interface described in `agent/outer.py`.
	"""

	def __init__(self, volume, concentrations):
		self.global_time = 0
		self.run_for = 1
		self.volume = volume
		self.simulations = {}
		self.concentrations = concentrations

	def time(self):
		return self.global_time

	def add_simulation(self, agent_id, state):
		self.simulations[agent_id] = state

	def remove_simulation(self, agent_id):
		return self.simulations.pop(agent_id, {})

	def update_from_simulations(self, update, now):
		self.simulations.update(update)

		for agent_id, state in update.iteritems():
			for molecule, change in state['changes'].iteritems():
				self.concentrations[molecule] += change

	def run_incremental(self, run_until):
		time.sleep(2)
		self.global_time = run_until

	def get_molecule_ids(self):
		return self.concentrations.keys()

	def get_concentrations(self, now):
		state = {}
		for agent_id, simulation in self.simulations.iteritems():
			if simulation['time'] <= now:
				state[agent_id] = self.concentrations

		return state

