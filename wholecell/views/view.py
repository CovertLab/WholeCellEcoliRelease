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


class ChromosomeMoleculeView(View):
	_stateID = 'Chromosome'

	def __init__(self, *args, **kwargs):
		super(ChromosomeMoleculeView, self).__init__(*args, **kwargs)


	# def _updateQuery(self):
	# 	self._state._container.

	
	def moleculeNew(self, moleculeName, location, **attributes):
		# if in partitioned region

		self._state._container.moleculeNew(moleculeName, location, 
			**attributes)


	def moleculeDel(self, molecule):
		# if in partitioned region

		self._state._container.moleculeDel(molecule)


	def moleculesDel(self, molecules):
		# if in partitioned region

		self._state._container.moleculesDel(molecules)


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
