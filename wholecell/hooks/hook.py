"""
Classes used to execute arbitrary code during critical parts of the simulation.
These are primarily used to support simulations that require behavior that 
cannot be modeled at this time.
"""

from __future__ import absolute_import, division, print_function


class SimulationHook(object):
	_name = 'SimulationHook'

	def __init__(self):
		pass


	@classmethod
	def name(cls):
		return cls._name


	def initialize(self, sim, sim_data):
		pass


	def postCalcInitialConditions(self, sim):
		pass


	def preEvolveState(self, sim):
		pass


	def postEvolveState(self, sim):
		pass


	def finalize(self, sim):
		pass

