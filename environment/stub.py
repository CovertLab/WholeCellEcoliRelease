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

	def time(self):
		return self.local_time

	def initialize_local_environment(self):
		pass

	def set_local_environment(self, molecule_ids, concentrations):
		self.molecule_ids = molecule_ids
		self.concentrations = concentrations

		if not self.local_set:
			self.environment_change = {}
			for molecule in molecule_ids:
				self.environment_change[molecule] = 0
			self.local_set = True

	def run_incremental(self, run_until):
		time.sleep(1)
		for molecule in self.molecule_ids:
			self.environment_change[molecule] += random.randint(1, 6)
		self.local_time = run_until

	def get_environment_change(self):
		return self.environment_change

	def finalize(self):
		pass
