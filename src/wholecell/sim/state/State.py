#!/usr/bin/env python

"""
State

State variable base class. Defines the interface states expose to the simulation and processes.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

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

		self.partitions = []

		for prop in propVals.keys():
			setattr(self, prop, propVals[prop])

	# Construct state-process graph, calculate constants
	def initialize(self, sim, kb):
		self.randStream = sim.randStream

	# Allocate memory
	def allocate(self):
		for partition in self.partitions:
			partition.allocate()

	# Calculate initial conditions
	def calcInitialConditions(self):
		return


	# -- Partitioning --

	def addPartition(self, process):
		try:
			partition = self.partitionClass(self, process)

		except:
			if self.partitionClass is None:
				raise Exception(
					'Called {0}.addPartition, but {0} has no partition class.'.format(type(self))
					)

			else:
				raise
		
		self.partitions.append(partition)

		return partition

	def prepartition(self):
		pass

	def partition(self):
		pass

	def merge(self):
		pass


	# -- Calculations --
	# Calculate (and cache) any dependent properties
	def calculate(self):
		return

	# -- Options and parameters

	def getOptions(self):
		val = {}
		if self.meta.has_key("options"):
			for opt in self.meta["options"]:
				val[opt] = getattr(self, opt)
		return val

	def setOptions(self, val):
		keys = val.keys()
		if not self.meta.has_key("options") or not all(set(keys).issubset(set(self.meta["options"]))):
			raise Exception, "Invalid option"

		for key in keys:
			setattr(self, key, val[key])

	def getParameters(self):
		val = {}
		if self.meta.has_key("parameters"):
			for param in self.meta["parameters"]:
				val[param] = getattr(self, param)
		return val

	def setParameters(self, val):
		keys = val.keys()
		if not self.meta.has_key("parameters") or not all(set(keys).issubset(set(self.meta["parameters"]))):
			raise Exception, "Invalid parameter"

		for key in keys:
			setattr(self, key, val[key])

	def getDynamics(self):
		val = {}
		for prop in self.meta["dynamics"]:
			val[prop] = getattr(self, prop)
		return val

	def setDynamics(self, val):
		keys = val.keys()
		if not self.meta.has_key("dynamics") or not all(set(keys).issubset(set(self.meta["dynamics"]))):
			raise Exception, "Invalid dynamics"

		for key in keys:
			setattr(self, key, val[key])
