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

	# Constructor
	def __init__(self, propVals = {}):
		# Metadata: id, name, list of dynamic properties, units
		if not hasattr(self, "meta"):
			self.meta = {}

		self._nProcesses = None
		self._views = None

		for prop in propVals.keys():
			setattr(self, prop, propVals[prop])


	# Construct state-process graph, calculate constants
	def initialize(self, sim, kb):
		self.randStream = sim.randStream

		self._nProcesses = len(sim.processes)
		self._views = []


	# Allocate memory
	def allocate(self):
		pass


	# Calculate initial conditions
	def calcInitialConditions(self):
		return


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

