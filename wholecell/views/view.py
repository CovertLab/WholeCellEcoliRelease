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

	def __init__(self, sim, process, query, evaluateQuery = True): # weight, priority, coupling key
		self._state = sim.states[self._stateID]
		self._state._queries.append(self)
		# self._process = process
		# self._processIndex = process._simulationIndex # TODO

		self._query = query # an immutable, hashable, composed of basic types
		self._evaluateQuery = evaluateQuery # whether to evaluate the query, which is not needed if the view is output-only


class BulkMoleculesView(View):
	_stateID = 'BulkMolecules'

	def __init__(self, *args, **kwargs):
		super(BulkMoleculesView, self).__init__(*args, **kwargs)

		# State references
		self._containerIndexes = self._state._container._namesToIndexes(self._query)

		# Memory allocation
		size = len(self._query)

		self._countsQuery = np.zeros(size, np.uint64) # TODO: revise names (countsTotal, countsRequested, countsAllocated?  drop "counts"?)
		self._countsRequest = np.zeros(size, np.uint64)
		self._counts = np.zeros(size, np.uint64)


	# Interface to State

	# TODO

	# Interface to Process

	# Query
	def countsQuery(self):
		return self._countsQuery.copy()

	# Request
	def countsRequestIs(self, values):
		# assert (self.countsQuery() >= values).all(), 'Request exceeds total values.'
		self._countsRequest[:] = values

	# Allocation
	def counts(self):
		return self._counts[:]


	def countsIs(self, values):
		# assert (values >= 0).all()
		self._counts[:] = values


	def countsInc(self, values):
		self._counts += values


	def countsDec(self, values):
		self._counts -= values

