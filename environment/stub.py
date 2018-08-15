from __future__ import absolute_import, division, print_function

import time
import random

class SimulationStub(object):

	"""
	Provide a stub for the simulation.

	This class is not meant to be functional but rather to implement in the simplest way 
	the interface the agents expect when interacting with the simulation.

	Full interface documented in `environment/inner.py`.
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

		if not self.local_set:
			self.environment_change = {}
			for molecule in self.concentrations.keys():
				self.environment_change[molecule] = 0
			self.local_set = True

	def run_incremental(self, run_until):
		time.sleep(1)
		for molecule in self.concentrations.keys():
			self.environment_change[molecule] += random.randint(1, 6)
		self.local_time = run_until

	def get_environment_change(self):
		return self.environment_change

	def finalize(self):
		pass

class EnvironmentStub(object):

	"""
	Provide a stub for the environmental context.

	Like the stub above for simulations, this class is meant to demonstrate the interface
	any environmental simulation must have to fulfill the conditions of being a coordinating
	external context for integrating the changes of each simulation.

	Full interface described in `environment/outer.py`.
	"""

	def __init__(self, volume, concentrations):
		self._time = 0
		self.run_for = 1
		self.volume = volume
		self.simulations = {}
		self.concentrations = concentrations

	def time(self):
		return self._time

	def add_simulation(self, id):
		state = {}
		self.simulations[id] = state

	def remove_simulation(self, id):
		return self.simulations.pop(id, {})

	def update_concentrations(self, all_changes):
		self._time += self.run_for
		for id, changes in all_changes.iteritems():
			for molecule, change in changes.iteritems():
				self.concentrations[molecule] += change

	def run_until(self):
		until = {}
		for id in self.simulations.keys():
			until[id] = self.time() + self.run_for

		return until

	def molecule_ids(self):
		return self.concentrations.keys()

	def get_concentrations(self):
		concentrations = {}
		for id in self.simulations.keys():
			concentrations[id] = self.concentrations
		
		return concentrations
