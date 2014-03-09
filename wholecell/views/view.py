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

		self._query = query # an immutable, hashable, composed of basic types

		self._totalCount = np.zeros(self._dataSize(), np.uint64) # number of objects that satisfy the query
		self._requestedCount = np.zeros_like(self._totalCount) # number of objects requested


	def _dataSize(self):
		return 1

	# Interface to State

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

	def _counts(self):
		return self._state._countsAllocated[self._containerIndexes, self._processIndex].copy()


	def _countsIs(self, values):
		self._state._countsAllocated[self._containerIndexes, self._processIndex] = values


	def _countsInc(self, values):
		self._state._countsAllocated[self._containerIndexes, self._processIndex] += values


	def _countsDec(self, values):
		self._state._countsAllocated[self._containerIndexes, self._processIndex] -= values


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

	# TODO
