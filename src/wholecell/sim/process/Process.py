#!/usr/bin/env python

"""
Process

Process submodel base class. Defines interface that processes expose to the simulation and to the states.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import inspect

class Process(object):
	""" Process """

	meta = {}

	# Constructor
	def __init__(self):
		# Constants
		self.timeStepSec = None

		# Simulation random stream
		self.randStream = None

	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		self.timeStepSec = sim.timeStepSec
		self.randStream = sim.randStream

	# Calculate submodel contribution to temporal evolution of cell
	def evolveState(self):
		return

	# -- Get, set options, parameters
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
