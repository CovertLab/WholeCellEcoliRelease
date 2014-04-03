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


class BulkMoleculesViewBase(View):
	_stateID = 'BulkMolecules'

	# TODO: implement update query method

	def _counts(self):
		return self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex].copy()


	def _countsIs(self, values):
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] = values


	def _countsInc(self, values):
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] += values


	def _countsDec(self, values):
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] -= values


class BulkMoleculesView(BulkMoleculesViewBase):
	def __init__(self, *args, **kwargs):
		super(BulkMoleculesView, self).__init__(*args, **kwargs)

		# State references
		self._containerIndexes = self._state._container._namesToIndexes(self._query)


	def _dataSize(self):
		return len(self._query)


	def counts(self):
		return self._counts()


	def countsIs(self, values):
		self._countsIs(values)


	def countsInc(self, values):
		self._countsInc(values)


	def countsDec(self, values):
		self._countsDec(values)


class BulkMoleculeView(BulkMoleculesViewBase):
	def __init__(self, *args, **kwargs):
		super(BulkMoleculeView, self).__init__(*args, **kwargs)

		# State references
		self._containerIndexes = self._state._container._namesToIndexes((self._query,))


	def count(self):
		return self._counts()


	def countIs(self, value):
		self._countsIs(value)


	def countInc(self, value):
		self._countsInc(value)


	def countDec(self, value):
		self._countsDec(value)


class UniqueMoleculesView(View):
	_stateID = 'UniqueMolecules'

	def __init__(self, *args, **kwargs):
		super(UniqueMoleculesView, self).__init__(*args, **kwargs)

		self._queryResult = None


	def _updateQuery(self):
		# TODO: generalize this logic (both here and in the state)

		self._queryResult = self._state.container.objectsWithName(
			self._query[0],
			**self._query[1]
			)

		self._totalIs(len(self._queryResult))


	def molecules(self):
		return self._state.container.objectsWithName(
			self._query[0],
			_partitionedProcess = ('==', self._processIndex + 1),
			**self._query[1]
			)

	# NOTE: these accessors do not enforce any sort of consistency between the query
	# and the objects created/deleted.  As such it may make more sense for these
	# to be process methods, not view methods. - JM
	def moleculeDel(self, molecule):
		self._state.container.objectDel(molecule)


	def moleculesDel(self, molecules):
		self._state.container.objectsDel(molecules)

	
	def moleculeNew(self, moleculeName, **attributes):
		self._state.container.objectNew(
			moleculeName,
			_partitionedProcess = ('==', self._processIndex + 1),
			**attributes
			)


	def moleculesNew(self, moleculeName, nMolecules, **attributes):
		self._state.container.objectsNew(
			moleculeName,
			nMolecules,
			_partitionedProcess = self._processIndex + 1,
			**attributes
			)


# TODO: views for forks
# TODO: determine how to assign regions to processes
# _array-sized matric with indexes on specific bases
# sparse representation (strand, start, stop)
# partition regions as unique objects?
# algo to merge overlapping regions (sort and perform comparisons on adjacent entries only should be O(n) + O(sort)) *use heapsort
# TODO: figure out interface methods, and how to constrain the regions that are operated in
# TODO: collect types of chrom views and organize in a base class
# TODO: views/operations for generic extents (i.e. 50-wide regions for binding)
# TODO: placeholder partitioning algo


class ChromosomeForksView(View):
	_stateID = 'Chromosome'

	# TODO: organize this class logically
	# TODO: incrementally add more methods
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

	# WARNING: for now, all queried regions can be operated on! no partitions!

	def parentRegions(self):
		return self._queryResult[0]


	def forksInRegion(self, regionParent):
		return self._state.container.forksInRegion(regionParent)


	def moleculeBoundOnFork(self, fork):
		# TODO: raise a specific exception
		# TODO: instead of constantly checking this, give each molecule/fork partitioned an attribute that can be checked
		assert self._state.container.forkInRegionSet(fork, self._queryResult[0])

		return self._state.container.moleculeBoundOnFork(fork)


	def maximumExtentPastFork(self, fork, maxExtentForward):
		assert self._state.container.forkInRegionSet(fork, self._queryResult[0])

		return self._state.container.maximumExtentPastFork(fork, maxExtentForward,
			self._queryResult[0])


	def moleculesBoundPastFork(self, fork, extentForward):
		molecules = self._state.container.moleculesBoundPastFork(fork, extentForward)

		assert all(
			self._state.container.moleculeInRegionSet(molecule, self._queryResult[0])
			for molecule in molecules
			)

		return molecules


	def moleculeLocationIsUnbound(self, molecule):
		# TODO: assert molecule can be accessed
		self._state.container.moleculeLocationIsUnbound(molecule)


	def forkExtend(self, fork, extent):
		assert self._state.container.forkInRegionSet(fork, self._queryResult[0])
		# TODO: assert 'extent' in region set

		newPosition = self._state.container.forkExtend(fork, extent)

		# TODO: update relevant regions in parent/children region sets


	def moleculeLocationIsFork(self, molecule, fork, extentForward, extentReverse):
		# TODO: assert region is accessible
		# TODO: assert molecule ownership by process

		self._state.container.moleculeLocationIsFork(molecule, fork,
			extentForward, extentReverse)

