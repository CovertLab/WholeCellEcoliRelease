#!/usr/bin/env python

"""
BulkMolecules.py

State which represents for a class of molecules the bulk copy numbers.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import re

import numpy as np

import wholecell.states.state
import wholecell.views.view
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

from wholecell.utils import units

from wholecell.utils.constants import REQUEST_PRIORITY_DEFAULT

ASSERT_POSITIVE_COUNTS = True

class BulkMolecules(wholecell.states.state.State):
	_name = 'BulkMolecules'

	def __init__(self, *args, **kwargs):
		self.container = None

		self._moleculeMass = None

		self._moleculeIDs = None
		# self._compartmentIDs = None

		# self._nCompartments = None

		self._countsRequested = None
		self._countsAllocatedInitial = None
		self._countsAllocatedFinal = None
		self._countsUnallocated = None

		super(BulkMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(BulkMolecules, self).initialize(sim, kb)

		# Load constants
		self._moleculeIDs = kb.bulkMolecules['moleculeId']
		# self._compartmentIDs = kb.compartments['compartmentAbbreviation']
		# self._nCompartments = kb.nCompartments

		self._moleculeMass = kb.bulkMolecules['mass'].asNumber(units.fg / units.mol) / kb.nAvogadro.asNumber(1 / units.mol)

		self._submassNameToIndex = kb.submassNameToIndex

		# self._compIndexes = {
		# 	compartmentKey:(kb.bulkMolecules['compartment'] == compartmentKey)
		# 	for compartmentKey in kb.compartments['compartmentAbbreviation']
		# 	}

		# Create the container for molecule counts
		self.container = BulkObjectsContainer(self._moleculeIDs)

		# Set up vector of process priorities
		self._processPriorities = np.empty(self._nProcesses, np.int64)
		self._processPriorities.fill(REQUEST_PRIORITY_DEFAULT)


	def processRequestPriorityIs(self, processIndex, priorityLevel):
		self._processPriorities[processIndex] = priorityLevel


	def allocate(self):
		super(BulkMolecules, self).allocate() # Allocates partitions

		nMolecules = self.container._counts.size
		dtype = self.container._counts.dtype

		# Arrays for tracking values related to partitioning
		self._countsRequested = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsAllocatedInitial = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsAllocatedFinal = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsUnallocated = np.zeros(nMolecules, dtype)


	def partition(self):
		if self._nProcesses == 0:
			self._countsUnallocated = self.container._counts
			return

		# Calculate and store requests
		self._countsRequested[:] = 0

		for view in self._views:
			self._countsRequested[view._containerIndexes, view._processIndex] += view._request()

		if ASSERT_POSITIVE_COUNTS:
			assert (self._countsRequested >= 0).all()

		# Calculate partition

		calculatePartition(self._processPriorities, self._countsRequested, self.container._counts, self._countsAllocatedInitial)

		if ASSERT_POSITIVE_COUNTS:
			assert (self._countsAllocatedInitial >= 0).all()

		# Record unpartitioned counts for later merging
		self._countsUnallocated = self.container._counts - np.sum(self._countsAllocatedInitial, axis = -1)

		if ASSERT_POSITIVE_COUNTS:
			assert (self._countsUnallocated >= 0).all()

		self._countsAllocatedFinal[:] = self._countsAllocatedInitial


	def calculatePreEvolveStateMass(self):
		# Compute masses of partitioned molecules

		if self.timeStep() == 0:
			self._countsUnallocated = self.container._counts

		self._masses[self._preEvolveStateMassIndex, ...] = np.dot(
			np.hstack([self._countsAllocatedInitial, self._countsUnallocated[:, np.newaxis]]).T,
			self._moleculeMass
			)


	def merge(self):
		if ASSERT_POSITIVE_COUNTS:
			assert (self._countsAllocatedFinal >= 0).all()

		self.container.countsIs(
			self._countsUnallocated + self._countsAllocatedFinal.sum(axis = -1)
			)


	def calculatePostEvolveStateMass(self):
		# Compute masses of partitioned molecules

		if self.timeStep() == 0:
			self._countsUnallocated = self.container._counts

		self._masses[self._postEvolveStateMassIndex, ...] = np.dot(
			np.hstack([self._countsAllocatedFinal, self._countsUnallocated[:, np.newaxis]]).T,
			self._moleculeMass
			)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			moleculeIDs = self._moleculeIDs
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			counts = self.container._counts,
			# TODO: better logging for requests i.e. not GBs of data
			# countsRequested = self._countsRequested,
			# countsAllocatedInitial = self._countsAllocatedInitial,
			# countsAllocatedFinal = self._countsAllocatedFinal,
			# countsUnallocated = self._countsUnallocated,
			)


	def tableLoad(self, tableReader, tableIndex):
		# entry = h5file.get_node('/', self._name)[timePoint]

		# self.container.countsIs(entry['counts'])

		# if self._nProcesses:
		# 	self._countsRequested[:] = entry['countsRequested']
		# 	self._countsAllocatedInitial[:] = entry['countsAllocatedInitial']
		# 	self._countsAllocatedFinal[:] = entry['countsAllocatedFinal']
		# 	self._countsUnallocated[:] = entry['countsUnallocated']

		raise NotImplementedError()


def calculatePartition(processPriorities, countsRequested, counts, countsPartitioned):
	# TODO: reduce the arrays to elements where counts != 0

	counts = counts.copy()

	priorityLevels = np.sort(np.unique(processPriorities))[::-1]

	for priorityLevel in priorityLevels:
		processHasPriority = (priorityLevel == processPriorities)

		requests = countsRequested[:, processHasPriority]

		totalRequests = requests.sum(axis = 1)
		totalRequestIsNonzero = (totalRequests > 0)

		fractionalRequests = np.zeros(requests.shape, np.float64)
		fractionalRequests[totalRequestIsNonzero, :] = (
			requests[totalRequestIsNonzero, :]
			/ totalRequests[totalRequestIsNonzero, np.newaxis]
			)

		allocations = np.fmin(
			requests,
			counts[:, np.newaxis] * fractionalRequests
			).astype(np.int64)

		countsPartitioned[:, processHasPriority] = allocations

		counts -= allocations.sum(axis = 1)

class BulkMoleculesViewBase(wholecell.views.view.View):
	_stateID = 'BulkMolecules'

	def _updateQuery(self):
		self._totalIs(self._state.container._counts[self._containerIndexes])


	def _counts(self):
		return self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex].copy()


	def _countsIs(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] = values


	def _countsInc(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] += values


	def _countsDec(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] -= values


class BulkMoleculesView(BulkMoleculesViewBase):
	def __init__(self, *args, **kwargs):
		super(BulkMoleculesView, self).__init__(*args, **kwargs)

		# State references
		assert len(set(self._query)) == len(self._query), "Bulk molecules views cannot contain duplicate entries"
		self._containerIndexes = self._state.container._namesToIndexes(self._query)


	def _dataSize(self):
		return len(self._query)


	def counts(self):
		return self._counts()


	def mass(self):
		return self._mass()


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
		self._containerIndexes = self._state.container._namesToIndexes((self._query,))


	def count(self):
		return self._counts()[0]


	def mass(self):
		return self._mass()


	def countIs(self, value):
		self._countsIs(value)


	def countInc(self, value):
		self._countsInc(value)


	def countDec(self, value):
		self._countsDec(value)

