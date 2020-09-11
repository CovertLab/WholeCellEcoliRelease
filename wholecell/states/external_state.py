"""
External State

State variable base class. Defines the interface states expose to the simulation.

"""

from __future__ import absolute_import, division, print_function


class ExternalState(object):
	""" External State """

	_name = 'ExternalState'

	# Constructor
	def __init__(self, *args, **kwargs):
		# Reference to sim
		self._sim = None

		# References to views
		self._views = []

		# Random number stream
		self.randomState = None

		self.seed = None


	# Construct state-process graph, calculate constants
	def initialize(self, sim, sim_data, timeline):
		self._sim = sim


	# Allocate memory
	def allocate(self):
		pass


	# Views
	def viewAdd(self, view):
		self._views.append(view)


	# Saving

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		pass


	# Basic accessors

	def time(self):
		return self._sim.time()

	def simulationStep(self):
		return self._sim.simulationStep()

	@classmethod
	def name(cls):
		return cls._name

