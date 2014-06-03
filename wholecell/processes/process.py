#!/usr/bin/env python

"""
Process

Process submodel base class. Defines interface that processes expose to the simulation and to the states.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import wholecell.views.view

import wholecell.states.bulk_molecules
import wholecell.states.unique_molecules

class Process(object):
	""" Process """

	_name = None

	# Constructor
	def __init__(self):
		# Constants
		self.timeStepSec = None
		self._processIndex = None

		# Simulation reference (used to access time)
		self._sim = None

		# Simulation random stream
		self.randomState = None # set by Simulation
		self.seed = None # set by Simulation

		# References to state
		self._states = None


	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		self._sim = sim

		self.timeStepSec = sim.timeStepSec
		self._processIndex = sim.processes.keys().index(self._name)

		self._states = sim.states


	# Construct views
	def bulkMoleculesView(self, moleculeIDs):
		return wholecell.states.bulk_molecules.BulkMoleculesView(
			self._states['BulkMolecules'], self, moleculeIDs)


	def bulkMoleculeView(self, moleculeIDs):
		return wholecell.states.bulk_molecules.BulkMoleculeView(
			self._states['BulkMolecules'], self, moleculeIDs)


	def uniqueMoleculesView(self, moleculeName, **attributes):
		return wholecell.states.unique_molecules.UniqueMoleculesView(
			self._states['UniqueMolecules'], self, (moleculeName, attributes))


	def chromosomeForksView(self, extentForward, extentReverse, includeMoleculesOnEnds):
		return wholecell.views.view.ChromosomeForksView(
			self._states['Chromosome'], self,
			(extentForward, extentReverse, includeMoleculesOnEnds))


	def chromosomeMoleculesView(self, moleculeName, extentForward, extentReverse, includeMoleculesOnEnds):
		return wholecell.views.view.ChromosomeMoleculesView(
			self._states['Chromosome'], self,
			(moleculeName, extentForward, extentReverse, includeMoleculesOnEnds))


	# Calculate requests for a single time step
	def calculateRequest(self):
		# Implemented by subclass
		pass


	# Calculate submodel contribution to temporal evolution of cell
	def evolveState(self):
		# Implemented by subclass
		pass


	# Basic accessors

	def time(self):
		return self._sim.time()


	@classmethod
	def name(cls):
		return cls._name
