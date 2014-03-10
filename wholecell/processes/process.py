#!/usr/bin/env python

"""
Process

Process submodel base class. Defines interface that processes expose to the simulation and to the states.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

class Process(object):
	""" Process """

	# Constructor
	def __init__(self):
		if not hasattr(self, "meta"):
			self.meta = {}

		# Constants
		self.timeStepSec = None
		self._processIndex = None

		# Simulation random stream
		self.randStream = None

		# References to state
		self._bulkMolecules = None


	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		self.timeStepSec = sim.timeStepSec
		self._processIndex = sim.processes.keys().index(self.meta['id'])

		self.randStream = sim.randStream

		self._bulkMolecules = sim.states['BulkMolecules']
		self._uniqueMolecules = sim.states['UniqueMolecules']


	# Construct views
	def bulkMoleculesView(self, moleculeIDs):
		import wholecell.views.view

		return wholecell.views.view.BulkMoleculesView(self._bulkMolecules, 
			self, moleculeIDs)


	def bulkMoleculeView(self, moleculeIDs):
		import wholecell.views.view

		return wholecell.views.view.BulkMoleculeView(self._bulkMolecules, 
			self, moleculeIDs)


	def uniqueMoleculesView(self, moleculeName, **attributes):
		import wholecell.views.view

		return wholecell.views.view.UniqueMoleculesView(self._uniqueMolecules,
			self, (moleculeName, attributes))


	# Calculate requests for a single time step
	def calculateRequest(self):
		pass


	# Calculate submodel contribution to temporal evolution of cell
	def evolveState(self):
		return


	# Partition requests
	def requestBulkMolecules(self):
		pass


	# -- Get, set options, parameters
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
