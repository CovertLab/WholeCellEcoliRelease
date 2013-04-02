#!/usr/bin/env python

"""
State

State variable base class. Defines the interface states expose to the simulation and processes.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import inspect

class State(object):
	""" State """

	# Metadata: id, name, list of dynamic properties, units
	meta = {}

	# Constructor
	def __init__(self, propVals = {}):
		# -- Partitioning --
		# Used by parent
		self.partitions = []

		# Used by children
		self.parentState = None
		self.parentProcess = None

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
		partition = self.constructPartition(process)
		self.partitions.append(partition)
		return partition

	def constructPartition(self, process):
		propVals = {}

		# TODO: May need to modify this depending on how concrete states are implemented (e.g. dependent variables)
		for prop, data in inspect.getmembers(self):
			if not callable(data):
				propVals[name] = data

		propVals["meta"]["id"] = "%s_%s" % (propVals["meta"]["id"], process["meta"]["id"])
		propVals["meta"]["name"] = "%s - %s" % (propVals["meta"]["name"], process["meta"]["name"])

		propVals["partitions"] = []
		propVals["parentState"] = self
		propVals["parentProcess"] = process

		return type(self)(propVals)

	def prepartition(self):
		return

	# Partition state among processes
	def partition(self):
		for partition in self.partitions:
			for dynamic in self.meta["dynamics"]:
				setattr(partition, dynamic, getattr(self, dynamic))

	# Merge sub-states partitioned to processes
	def merge(self):
		for dynamic in self.meta["dynamics"]:
			oldVal = getattr(self, dynamic)
			newVal = oldVal
			nNewVal = 0

			for partition in self.partitions:
				if oldVal != getattr(partition, dynamic):
					newVal = getattr(partition, dynamic)
					nNewVal += 1

			if nNewVal > 1:
				raise Exception, "Multiple processes cannot simultaneously edit state properties"

			setattr(self, dynamic, newVal)


	# -- Calculations --
	# Calculate (and cache) any dependent properties
	def calculate(self):
		return

	# -- Options and parameters

	def getOptions(self):
		val = {}
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

	def setDynamics(self, val):
		keys = val.keys()
		if not self.meta.has_key("dynamics") or not all(set(keys).issubset(set(self.meta["dynamics"]))):
			raise Exception, "Invalid dynamics"

		for key in keys:
			setattr(self, key, val[key])