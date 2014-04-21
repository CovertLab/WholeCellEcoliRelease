'''
view.py

A View is an interface between a State and a Process.  It provides support for
querying the State, requesting part of the State, providing the allocation to 
the Process, and returning the updated values to the State for merging. The 
implementation is largely specific to the State.
'''

from __future__ import division

import numpy as np

class View(object):
	_stateID = None

	def __init__(self, state, process, query): # weight, priority, coupling id, option to not evaluate the query
		self._state = state
		self._state.viewAdd(self)
		self._processId = process.name()
		self._processIndex = process._processIndex

		self._query = query

		self._totalCount = np.zeros(self._dataSize(), np.int64) # number of objects that satisfy the query
		self._requestedCount = np.zeros_like(self._totalCount) # number of objects requested


	def _dataSize(self):
		return 1

	# Interface to State

	def _updateQuery(self):
		pass


	def _totalIs(self, value):
		self._totalCount[:] = value

		# Clear out request in preparation
		self._requestedCount[:] = 0


	def _request(self):
		return self._requestedCount # NOTE: this is not a copy - be careful!


	# Interface to Process

	# Query

	def total(self):
		return self._totalCount.copy()

	# Request

	def requestIs(self, value):
		assert (value <= self._totalCount).all(), 'Requested more than exists'
		self._requestedCount[:] = value


	def requestAll(self):
		self._requestedCount[:] = self._totalCount

# TODO: views/operations for generic extents (i.e. 50-wide regions for binding)

from wholecell.containers.chromosome_container import _ChromosomeRegionSet # TODO: don't import a private class

class ChromosomeForksView(View):
	_stateID = 'Chromosome'

	# TODO: organize this class logically
	# TODO: move most methods to a base ChromosomeView class

	def _updateQuery(self):
		# Query structure:
		# extentForward
		# extentReverse
		# includeMoleculesOnEnds

		self._queryResult = self._state.container.regionsNearForks(
			*self._query
			)

		self._totalIs(len(self._queryResult[0]))


	def requestedRegions(self):
		for regionSet in self._queryResult:
			for region in regionSet:
				yield region


	def allocateRegions(self, regions):
		# NOTE: this code is really bad
		queryRegionsParent = list(self._queryResult[0])
		queryRegionsChildA = list(self._queryResult[1])
		queryRegionsChildB = list(self._queryResult[2])

		allocatedRegionsParent = []
		allocatedRegionsChildA = []
		allocatedRegionsChildB = []

		for region in regions:
			if region in queryRegionsParent:
				group = allocatedRegionsParent

			elif region in queryRegionsChildA:
				group = allocatedRegionsChildA

			elif region in queryRegionsChildB:
				group = allocatedRegionsChildB

			group.append(
				(region._strand, region._start, region._stop)
				)

		self._allocatedRegions = (
			_ChromosomeRegionSet(allocatedRegionsParent),
			_ChromosomeRegionSet(allocatedRegionsChildA),
			_ChromosomeRegionSet(allocatedRegionsChildB),
			)

		self._allocatedMolecules = (
			self._state.container.moleculesInRegionSet(self._allocatedRegions[0])
			| self._state.container.moleculesInRegionSet(self._allocatedRegions[1])
			| self._state.container.moleculesInRegionSet(self._allocatedRegions[2])
			)


	def parentRegions(self):
		return self._allocatedRegions[0]


	def forksInRegion(self, regionParent):
		return self._state.container.forksInRegion(regionParent)


	def moleculeBoundOnFork(self, fork):
		# TODO: raise a specific exception
		# TODO: instead of constantly checking this, give each molecule/fork partitioned an attribute that can be checked
		assert self._state.container.forkInRegionSet(fork, self._allocatedRegions[0])

		return self._state.container.moleculeBoundOnFork(fork)


	def maximumExtentPastFork(self, fork, maxExtentForward):
		assert self._state.container.forkInRegionSet(fork, self._allocatedRegions[0])

		return self._state.container.maximumExtentPastFork(fork, maxExtentForward,
			self._allocatedRegions[0])


	def moleculesBoundPastFork(self, fork, extentForward):
		molecules = self._state.container.moleculesBoundPastFork(fork, extentForward)

		assert all(
			self._state.container.moleculeInRegionSet(molecule, self._allocatedRegions[0])
			for molecule in molecules
			)

		return molecules


	def moleculeLocationIsUnbound(self, molecule):
		assert molecule in self._allocatedMolecules

		self._state.container.moleculeLocationIsUnbound(molecule)


	def forkExtend(self, fork, extent):
		assert self._state.container.forkInRegionSet(fork, self._allocatedRegions[0])
		# TODO: assert 'extent' in region set

		newPosition = self._state.container.forkExtend(fork, extent)

		# TODO: update relevant regions in parent/children region sets


	def moleculeLocationIsFork(self, molecule, fork, extentForward, extentReverse):
		# TODO: assert region is accessible
		# TODO: assert molecule ownership by process

		assert molecule in self._allocatedMolecules
		assert self._state.container.forkInRegionSet(fork, self._allocatedRegions[0])

		self._state.container.moleculeLocationIsFork(molecule, fork,
			extentForward, extentReverse)


class ChromosomeMoleculesView(View):
	_stateID = 'Chromosome'

	# TODO: reconcile with other chromosome view

	def _updateQuery(self):
		# Query structure:
		# molecule name
		# extentForward
		# extentReverse
		# includeMoleculesOnEnds

		self._queryResult = self._state.container.regionsNearMolecules(
			*self._query
			)

		self._totalIs(len(self._queryResult))


	def requestedRegions(self):
		for region in self._queryResult:
			yield region


	def allocateRegions(self, regions):
		allocatedRegions = []

		for region in regions:
			allocatedRegions.append(
				(region._strand, region._start, region._stop)
				)
			
		self._allocatedRegions = _ChromosomeRegionSet(allocatedRegions)

		self._allocatedMolecules = self._state.container.moleculesInRegionSet(self._allocatedRegions)


	def regions(self):
		return self._allocatedRegions


	def moleculeLocationIsUnbound(self, molecule):
		assert molecule in self._allocatedMolecules

		self._state.container.moleculeLocationIsUnbound(molecule)

