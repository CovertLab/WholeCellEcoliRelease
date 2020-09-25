"""
Internal State

State variable base class. Defines the interface states expose to the simulation and processes.

"""

from __future__ import absolute_import, division, print_function

import numpy as np

class InternalState(object):
	""" Internal State """

	_name = 'InternalState'

	# Constructor
	def __init__(self):
		# Constants
		self._nProcesses = None

		# Reference to sim
		self._sim = None

		# References to views
		self._views = []

		# Random number stream
		self.randomState = None  # Set at the first time step by Simulation._evolveState()
		self.seed = None

		# Mass arrays
		self._masses = None
		self._process_mass_diffs = None


	# Construct state-process graph, calculate constants
	def initialize(self, sim, sim_data):
		self._sim = sim

		self._nProcesses = len(sim.processes)

		self._masses = np.zeros(len(sim_data.submass_name_to_index), np.float64)
		self._compartment_masses = np.zeros(
			(len(sim_data.compartment_abbrev_to_index),
			 len(sim_data.submass_name_to_index)
			 ), np.float64)
		self._process_mass_diffs = np.zeros(
			(self._nProcesses, len(sim_data.submass_name_to_index)), np.float64)


	# Allocate memory
	def allocate(self):
		pass


	# Views
	def viewAdd(self, view):
		self._views.append(view)


	# Partitioning
	def updateQueries(self):
		for view in self._views:
			view._updateQuery()


	def partition(self, processes):
		pass


	def merge(self, processes):
		pass


	def calculateMass(self):
		raise NotImplementedError("Subclass must implement")


	# Mass calculations
	def mass(self):
		return self._masses

	# Mass calculations by compartment
	def compartment_mass(self):
		return self._compartment_masses

	def process_mass_diffs(self):
		return self._process_mass_diffs

	def reset_process_mass_diffs(self):
		self._process_mass_diffs.fill(0)


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

