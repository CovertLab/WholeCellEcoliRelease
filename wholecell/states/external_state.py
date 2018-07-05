#!/usr/bin/env python

"""
External State

State variable base class. Defines the interface states expose to the simulation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import numpy as np

class ExternalState(object):
	""" External State """

	_name = None

	# Constructor
	def __init__(self):
		# Reference to sim
		self._sim = None

		# Random number stream
		self.randomState = None

		self.seed = None


	# Construct state-process graph, calculate constants
	def initialize(self, sim, sim_data):
		self._sim = sim


	# Allocate memory
	def allocate(self):
		pass


	# Saving and loading

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		pass


	def tableLoad(self, tableReader, tableIndex):
		pass


	# Basic accessors

	def time(self):
		return self._sim.time()

	def simulationStep(self):
		return self._sim.simulationStep()

	@classmethod
	def name(cls):
		return cls._name

