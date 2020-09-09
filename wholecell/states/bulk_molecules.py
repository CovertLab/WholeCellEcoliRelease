#!/usr/bin/env python

"""
BulkMolecules.py

State which represents for a class of molecules the bulk copy numbers.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from six.moves import zip

import wholecell.states.internal_state
import wholecell.views.view
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils import units

from wholecell.utils.constants import REQUEST_PRIORITY_DEFAULT

ASSERT_POSITIVE_COUNTS = True


class NegativeCountsError(Exception):
	pass


class BulkMolecules(wholecell.states.internal_state.InternalState):
	_name = 'BulkMolecules'

	def __init__(self):
		super(BulkMolecules, self).__init__()

		self.container = None
		self._moleculeMass = None
		self._moleculeIDs = None
		self._countsRequested = None
		self._countsAllocatedInitial = None
		self._countsAllocatedFinal = None
		self._countsUnallocated = None

		self._processIDs = None
		self._processID_to_index = {}
		self._submass_name_to_index = None
		self._processPriorities = None
		self.division_mode = {}


	def initialize(self, sim, sim_data):
		super(BulkMolecules, self).initialize(sim, sim_data)

		self._processIDs = list(sim.processes.keys())
		self._processID_to_index = {
			id_: idx for idx, id_ in enumerate(self._processIDs)}

		# Load constants
		self._moleculeIDs = sim_data.internal_state.bulk_molecules.bulk_data['id']

		self._moleculeMass = sim_data.internal_state.bulk_molecules.bulk_data['mass'].asNumber(units.fg / units.mol) / sim_data.constants.n_avogadro.asNumber(1 / units.mol)

		self._submass_name_to_index = sim_data.submass_name_to_index

		# Create the container for molecule counts
		self.container = BulkObjectsContainer(self._moleculeIDs)

		# Set up vector of process priorities
		self._processPriorities = np.empty(self._nProcesses, np.int64)
		self._processPriorities.fill(REQUEST_PRIORITY_DEFAULT)

		# Set up ids for division into daughter cells
		self.division_mode = {}
		self.division_mode['binomial'] = sim_data.molecule_groups.bulk_molecules_binomial_division
		self.division_mode['equally'] = sim_data.molecule_groups.bulk_molecules_equal_division

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


	def partition(self, processes):
		# Reset allocated counts
		self._countsAllocatedInitial.fill(0)

		if len(processes) == 0:
			self._countsUnallocated = self.container._counts
			return

		# Get list of process indexes to be partitioned
		process_indexes = np.array([
			self._processID_to_index[process.__name__]
		    for process in processes])

		# Calculate and store requests
		self._countsRequested.fill(0)

		process_indexes_set = set(process_indexes)
		for view in self._views:
			if view._processIndex in process_indexes_set:
				self._countsRequested[view._containerIndexes, view._processIndex] += view._request()

		# Select columns to partition
		counts_requested = self._countsRequested[:, process_indexes]

		if ASSERT_POSITIVE_COUNTS and np.any(counts_requested < 0):
			raise NegativeCountsError(
				"Negative value(s) in self._countsRequested:\n"
				+ "\n".join(
					"{} in {} ({})".format(
						self._moleculeIDs[molIndex],
						self._processIDs[processIndex],
						self._countsRequested[molIndex, processIndex]
						)
					for molIndex, processIndex in zip(*np.where(self._countsRequested < 0))
					)
				)

		# Calculate partition
		if len(processes) > 1:
			self._countsAllocatedInitial[:, process_indexes] = calculatePartition(
				self._processPriorities[process_indexes],
				counts_requested,
				self.container._counts,
				self.randomState,
				)
		else:
			# No need to partition if there is only one process
			self._countsAllocatedInitial[:, process_indexes[0]] = np.fmin(
				counts_requested.flatten(), self.container._counts)

		# Select columns that have been partitioned
		counts_allocated_initial = self._countsAllocatedInitial[:, process_indexes]

		if ASSERT_POSITIVE_COUNTS and np.any(counts_allocated_initial < 0):
			raise NegativeCountsError(
					"Negative value(s) in self._countsAllocatedInitial:\n"
					+ "\n".join(
					"{} in {} ({})".format(
						self._moleculeIDs[molIndex],
						self._processIDs[processIndex],
						self._countsAllocatedInitial[molIndex, processIndex]
						)
					for molIndex, processIndex in zip(*np.where(self._countsAllocatedInitial < 0))
					)
				)

		# Record unpartitioned counts for later merging
		self._countsUnallocated = self.container._counts - np.sum(
			counts_allocated_initial, axis=-1)

		if ASSERT_POSITIVE_COUNTS and np.any(self._countsUnallocated < 0):
			raise NegativeCountsError(
					"Negative value(s) in self._countsUnallocated:\n"
					+ "\n".join(
					"{} ({})".format(
						self._moleculeIDs[molIndex],
						self._countsUnallocated[molIndex]
						)
					for molIndex in np.where(self._countsUnallocated < 0)[0]
					)
				)

		np.copyto(self._countsAllocatedFinal, self._countsAllocatedInitial)


	def merge(self, processes):
		if len(processes) == 0:
			return

		# Get list of process indexes to be merged
		process_indexes = np.array([
			self._processID_to_index[process.__name__]
			for process in processes])

		# Select columns to merge
		counts_allocated_final = self._countsAllocatedFinal[:, process_indexes]

		if ASSERT_POSITIVE_COUNTS and np.any(counts_allocated_final < 0):
			raise NegativeCountsError(
					"Negative value(s) in self._countsAllocatedFinal:\n"
					+ "\n".join(
					"{} in {} ({})".format(
						self._moleculeIDs[molIndex],
						self._processIDs[processIndex],
						self._countsAllocatedFinal[molIndex, processIndex]
						)
					for molIndex, processIndex in zip(*np.where(self._countsAllocatedFinal < 0))
					)
				)

		# Merge counts and calculate mass differences for each process
		if len(processes) > 1:
			self.container.countsIs(
				self._countsUnallocated + counts_allocated_final.sum(axis=-1)
				)
			self._process_mass_diffs[process_indexes, :] += np.dot(
				(self._countsAllocatedFinal - self._countsAllocatedInitial)[:, process_indexes].T,
				self._moleculeMass)
		else:
			# Use simpler calculations if there is only one process
			counts_allocated_final = counts_allocated_final.flatten()
			self.container.countsIs(
				self._countsUnallocated + counts_allocated_final
				)
			self._process_mass_diffs[process_indexes[0], :] += np.dot(
				counts_allocated_final - self._countsAllocatedInitial[:, process_indexes[0]],
				self._moleculeMass)


	def calculateMass(self):
		# Compute summed masses of all molecules
		if self.simulationStep() == 0:
			self._countsUnallocated = self.container._counts

		self._masses = np.dot(
			np.hstack([self._countsAllocatedFinal, self._countsUnallocated[:, np.newaxis]]).T,
			self._moleculeMass
			).sum(axis=0)


	def loadSnapshot(self, container):
		"""Load data from a snapshot `container`."""
		self.container.loadSnapshot(container)

	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		subcolumns = {
			'counts': 'objectNames'}

		tableWriter.writeAttributes(
			processNames = self._processIDs,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			counts = self.container._counts,
			atpAllocatedInitial = self._countsAllocatedInitial[self.container._objectNames.index("ATP[c]"), :],
			atpAllocatedFinal = self._countsAllocatedFinal[self.container._objectNames.index("ATP[c]"), :],
			atpRequested = self._countsRequested[self.container._objectNames.index("ATP[c]"), :],
			)


def calculatePartition(processPriorities, countsRequested, counts, random_state):
	partitioned_counts = np.zeros_like(countsRequested)
	counts = counts.copy()

	priorityLevels = np.sort(np.unique(processPriorities))[::-1]

	for priorityLevel in priorityLevels:
		processHasPriority = (priorityLevel == processPriorities)

		requests = countsRequested[:, processHasPriority].copy()

		totalRequests = requests.sum(axis=1)
		excess_request_mask = (totalRequests > counts)

		# Get fractional request for molecules that have excess request
		# compared to available counts
		fractional_requests = (
			requests[excess_request_mask, :] * counts[excess_request_mask, np.newaxis]
			/ totalRequests[excess_request_mask, np.newaxis]
			)

		# Distribute fractional counts to ensure full allocation of excess
		# request molecules
		remainders = fractional_requests % 1
		for idx, remainder in enumerate(remainders):
			count = int(np.round(remainder.sum()))
			fractional_requests[idx, :] += random_state.multinomial(count, remainder)
		requests[excess_request_mask, :] = fractional_requests

		allocations = requests.astype(np.int64)
		partitioned_counts[:, processHasPriority] = allocations
		counts -= allocations.sum(axis=1)

	return partitioned_counts


class BulkMoleculesViewBase(wholecell.views.view.View):
	_stateID = 'BulkMolecules'

	def __init__(self, *args, **kwargs):
		super(BulkMoleculesViewBase, self).__init__(*args, **kwargs)
		self._containerIndexes = None  # subclasses must set this

	def _updateQuery(self):
		self._totalIs(self._state.container._counts[self._containerIndexes])


	def _counts(self):
		return self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex].copy()


	def _countsIs(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] = values


	def _countsInc(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		values = np.asarray(values, dtype=self._containerIndexes.dtype)
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] += values


	def _countsDec(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		values = np.asarray(values, dtype=self._containerIndexes.dtype)
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] -= values

	# Request
	def requestIs(self, value):
		self._requestedCount[:] = value

	def requestAll(self):
		np.copyto(self._requestedCount, self._totalCount)


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


	def countIs(self, value):
		self._countsIs(value)


	def countInc(self, value):
		self._countsInc(value)


	def countDec(self, value):
		self._countsDec(value)

