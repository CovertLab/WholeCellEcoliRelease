"""
Classes used to execute arbitrary code during critical parts of the simulation.
These are primarily used to support simulations that require behavior that 
cannot be modeled at this time.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import division

import numpy as np

class SimulationHook(object):
	_name = None

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

