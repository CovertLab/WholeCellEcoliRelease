#!/usr/bin/env python

"""
State

State variable base class. Defines the interface states expose to the simulation and processes.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

from __future__ import division

class State(object):
	""" State """

	_name = None

	# Constructor
	def __init__(self):
		# Constants
		self._nProcesses = None

		# Reference to sim
		self._sim = None

		# References to views
		self._views = []

		# Random number stream
		self.randStream = None


	# Construct state-process graph, calculate constants
	def initialize(self, sim, kb):
		self._sim = sim

		self._nProcesses = len(sim.processes)


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


	def partition(self):
		pass


	def merge(self):
		pass


	# Saving and loading

	def pytablesCreate(self, h5file, expectedRows):
		pass

	def pytablesAppend(self, h5file):
		pass

	def pytablesLoad(self, h5file, timePoint):
		pass


	# Calculations

	# Calculate (and cache) any dependent properties
	def calculate(self):
		return


	# Basic accessors

	def timeStep(self):
		return self._sim.timeStep()

	@classmethod
	def name(cls):
		return cls._name

