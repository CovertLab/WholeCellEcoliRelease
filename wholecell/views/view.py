'''
view.py

A View is an interface between a State and a Process.  It provides support for
querying the State, requesting part of the State, providing the allocation to
the Process, and returning the updated values to the State for merging. The
implementation is largely specific to the State.
'''

from __future__ import absolute_import, division, print_function

import numpy as np

class TotalCountError(Exception):
	pass

class View(object):
	_stateID = 'View'

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

	def total_counts(self):
		return self._totalCount.copy()

	def total_count(self):
		if self._totalCount.size != 1:
			raise TotalCountError('total_count() can only be used for views that cover a single molecule type.')

		return self._totalCount[0]

	# TODO (ggsun): deprecated alias, should be deleted
	total = total_counts
