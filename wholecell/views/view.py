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

	def __init__(self, sim, process, query): # weight, priority, coupling id, option to not evaluate the query
		self._state = sim.states[self._stateID]
		sim.states[self._stateID]._queries.append(self)
		self._processIndex = process._simulationIndex

		self._query = query # an immutable, hashable, composed of basic types

		self._totalCount = np.zeros(self._dataSize(), np.uint64) # number of objects that satisfy the query
		self._requestedCount = np.zeros_like(self._totalCount) # number of objects requested


	def _dataSize(self):
		return 1

	# Interface to Process

	# Query

	def total(self):
		return self._totalCount.copy()

	# Request

	def requestIs(self, value):
		assert (value <= )
		self._requestedCount[:] = value


class BulkMoleculesView(View):
	_stateID = 'BulkMolecules'

	def __init__(self, *args, **kwargs):
		super(BulkMoleculesView, self).__init__(*args, **kwargs)

		# State references
		self._containerIndexes = self._state._container._namesToIndexes(self._query)

		# Memory allocation
		self._counts = np.zeros_like(self._totalNumber)


	def _dataSize(self):
		return len(self._query)


	# Interface to State

	# TODO

	# Interface to Process

	# TODO

