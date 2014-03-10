#!/usr/bin/env python

"""
State

State variable base class. Defines the interface states expose to the simulation and processes.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

from __future__ import division

from collections import OrderedDict

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
		pass


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

	# Options and parameters

	def getOptions(self):
		val = {}
		if self.meta.has_key("options"):
			for opt in self.meta["options"]:
				val[opt] = getattr(self, opt)
		return val

	def setOptions(self, val):
		keys = val.keys()
		if not self.meta.has_key("options") or not set(keys).issubset(set(self.meta["options"])):
			raise Exception, "Invalid option"

		for key in keys:
			setattr(self, key, val[key])

	def getDynamics(self):
		val = {}
		for prop in self.meta["dynamics"]:
			val[prop] = getattr(self, prop)
		return val