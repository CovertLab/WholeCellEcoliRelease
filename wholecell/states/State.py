#!/usr/bin/env python

"""
State

State variable base class. Defines the interface states expose to the simulation and processes.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

from collections import OrderedDict

import numpy

class State(object):
	""" State """

	partitions = None
	partitionClass = None

	# Constructor
	def __init__(self, propVals = {}):
		# Metadata: id, name, list of dynamic properties, units
		if not hasattr(self, "meta"):
			self.meta = {}

		self.partitions = OrderedDict()

		for prop in propVals.keys():
			setattr(self, prop, propVals[prop])

	# Construct state-process graph, calculate constants
	def initialize(self, sim, kb):
		self.randStream = sim.randStream

		if self.partitionClass is not None:
			for process in sim.processes.viewvalues():
				self.partitions[process] = self.partitionClass(self, process)

	# Allocate memory
	def allocate(self):
		for partition in self.partitions.viewvalues():
			partition.allocate()

	# Calculate initial conditions
	def calcInitialConditions(self):
		return


	# Partitioning
	def setPartition(self):
		pass

	def prepartition(self):
		pass

	def partition(self):
		pass

	def merge(self):
		pass


	# Saving and loading

	def pytablesCreate(self, h5file):
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