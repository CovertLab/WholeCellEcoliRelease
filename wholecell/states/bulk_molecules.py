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
import tables

import wholecell.states.state
import wholecell.views.view
from wholecell.containers.bulk_objects_container import BulkObjectsContainer


class BulkMolecules(wholecell.states.state.State):
	_name = 'BulkMolecules'

	def __init__(self, *args, **kwargs):
		self.container = None

		self._moleculeMass = None

		self._moleculeIDs = None
		self._compartmentIDs = None

		self._nCompartments = None

		self._isRequestAbsolute = None

		self._countsRequested = None
		self._countsAllocatedInitial = None
		self._countsAllocatedFinal = None
		self._countsUnallocated = None

		self._typeIdxs = None

		super(BulkMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(BulkMolecules, self).initialize(sim, kb)

		# Load constants
		self._moleculeIDs = moleculeIds(kb)
		self._compartmentIDs = kb.compartments['compartmentAbbreviation']
		self._nCompartments = kb.nCompartments

		self._moleculeMass = kb.bulkMolecules['mass'].to('fg / mol').magnitude / kb.nAvogadro.to('1 / mole').magnitude

		self._typeIdxs = {'metabolites'	:	kb.bulkMolecules['isMetabolite'],
							'rnas'		:	kb.bulkMolecules['isRnaMonomer'],
							'rrnas'		:	np.array([True if x in kb.rnaData["id"][kb.rnaData["isRRna"]] else False for x in kb.bulkMolecules["moleculeId"]]),
							'proteins'	:	kb.bulkMolecules['isProteinMonomer'],
							'water'		:	kb.bulkMolecules['isWater']}

		self._compIndexes = {
			compartmentKey:(kb.bulkMolecules['compartment'] == compartmentKey)
			for compartmentKey in kb.compartments['compartmentAbbreviation']
			}

		# Create the container for molecule counts
		self.container = bulkObjectsContainer(kb)
		
		# TODO: restore this behavior or replace it with something bettter

		self._isRequestAbsolute = np.zeros(self._nProcesses, np.bool)
		try:
			self._isRequestAbsolute[sim.processes['RnaDegradation']._processIndex] = True

		except KeyError:
			pass

	def allocate(self):
		super(BulkMolecules, self).allocate() # Allocates partitions

		nMolecules = self.container._counts.size
		dtype = self.container._counts.dtype

		# Arrays for tracking values related to partitioning
		self._countsRequested = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsAllocatedInitial = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsAllocatedFinal = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsUnallocated = np.zeros(nMolecules, dtype)


	def updateQueries(self):
		for view in self._views:
			view._totalIs(self.container._counts[view.containerIndexes])


	def partition(self):
		if self._nProcesses:
			# Calculate and store requests
			self._countsRequested[:] = 0

			for view in self._views:
				self._countsRequested[view.containerIndexes, view._processIndex] += view._request()

			calculatePartition(self._isRequestAbsolute, self._countsRequested, self.container._counts, self._countsAllocatedInitial)
			
			# Record unpartitioned counts for later merging
			self._countsUnallocated = self.container._counts - np.sum(self._countsAllocatedInitial, axis = -1)

			self._countsAllocatedFinal[:] = self._countsAllocatedInitial

			self._massAllocatedInitial = (
				self._countsAllocatedInitial *
				np.tile(self._moleculeMass.reshape(-1, 1), (1, self._nProcesses))
				)

		else:
			self._countsUnallocated = self.container._counts


	def merge(self):
		self.container.countsIs(
			self._countsUnallocated + self._countsAllocatedFinal.sum(axis = -1)
			)
		self._massAllocatedFinal = (
			self._countsAllocatedFinal *
			np.tile(self._moleculeMass.reshape(-1, 1), (1, self._nProcesses))
			)


	def mass(self):
		return np.dot(
			self._moleculeMass,
			self.container._counts
			)


	def massByType(self, typeKey):
		indexes = self._typeIdxs[typeKey]

		return np.dot(
			self._moleculeMass[indexes],
			self.container._counts[indexes]
			)


	def massByCompartment(self, compartment):
		indexes = self._compIndexes[compartment]

		return np.dot(
			self._moleculeMass[indexes],
			self.container._counts[indexes]
			)


	def pytablesCreate(self, h5file, expectedRows):
		countsShape = self.container._counts.shape
		partitionsShape = self._countsRequested.shape

		# Columns
		d = {
			"time": tables.Int64Col(),
			"counts":tables.UInt64Col(countsShape),
			"countsRequested":tables.UInt64Col(partitionsShape),
			"countsAllocatedInitial":tables.UInt64Col(partitionsShape),
			"countsAllocatedFinal":tables.UInt64Col(partitionsShape),
			"countsUnallocated":tables.UInt64Col(countsShape),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self._name,
			d,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)
	
		groupNames = h5file.create_group(h5file.root,
			'names', 'Molecule, compartment, and process names')

		h5file.create_array(groupNames, 'moleculeIDs', [str(s) for s in self._moleculeIDs]) # pytables doesn't support unicode
		h5file.create_array(groupNames, 'compartmentIDs', [str(s) for s in self._compartmentIDs])

		groupIdxs = h5file.create_group(h5file.root,
			'indexes', 'Indexes for various groups of molecules')

		for type_, indexes in self._typeIdxs.viewitems():
			h5file.create_array(groupIdxs, type_, indexes)


	def pytablesAppend(self, h5file):
		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.timeStep()
		entry['counts'] = self.container._counts
		entry['countsRequested'] = self._countsRequested
		entry['countsAllocatedInitial'] = self._countsAllocatedInitial
		entry['countsAllocatedFinal'] = self._countsAllocatedFinal
		entry['countsUnallocated'] = self._countsUnallocated
		
		entry.append()

		t.flush()


	def pytablesLoad(self, h5file, timePoint):
		entry = h5file.get_node('/', self._name)[timePoint]

		self.container.countsIs(entry['counts'])
		
		if self._nProcesses:
			self._countsRequested[:] = entry['countsRequested']
			self._countsAllocatedInitial[:] = entry['countsAllocatedInitial']
			self._countsAllocatedFinal[:] = entry['countsAllocatedFinal']
			self._countsUnallocated[:] = entry['countsUnallocated']


def calculatePartition(isRequestAbsolute, countsBulkRequested, countsBulk, countsBulkPartitioned):
	requestsAbsolute = np.sum(countsBulkRequested[..., isRequestAbsolute], axis = -1)
	requestsRelative = np.sum(countsBulkRequested[..., ~isRequestAbsolute], axis = -1)

	# TODO: Remove the warnings filter or move it elsewhere
	# there may also be a way to avoid these warnings by only evaluating 
	# division "sparsely", which should be faster anyway - JM
	oldSettings = np.seterr(invalid = 'ignore', divide = 'ignore') # Ignore divides-by-zero errors

	scaleAbsolute = np.fmax(0, # Restrict requests to at least 0% (fmax replaces nan's)
		np.minimum(1, # Restrict requests to at most 100% (absolute requests can do strange things)
			np.minimum(countsBulk, requestsAbsolute) / requestsAbsolute) # Divide requests amongst partitions proportionally
		)

	scaleRelative = np.fmax(0, # Restrict requests to at least 0% (fmax replaces nan's)
		np.maximum(0, countsBulk - requestsAbsolute) / requestsRelative # Divide remaining requests amongst partitions proportionally
		)

	scaleRelative[requestsRelative == 0] = 0 # nan handling?

	np.seterr(**oldSettings) # Restore error handling to the previous state

	# Compute allocations and assign counts to the partitions
	for iPartition in range(countsBulkPartitioned.shape[-1]):
		scale = scaleAbsolute if isRequestAbsolute[iPartition] else scaleRelative
		allocation = np.floor(countsBulkRequested[..., iPartition] * scale)
		countsBulkPartitioned[..., iPartition] = allocation


def moleculeIds(kb):
	return kb.bulkMolecules['moleculeId']


def bulkObjectsContainer(kb, dtype = np.int64):
	return BulkObjectsContainer(moleculeIds(kb), dtype)


class BulkMoleculesViewBase(wholecell.views.view.View):
	_stateID = 'BulkMolecules'

	def _counts(self):
		return self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex].copy()


	def _mass(self):
		return numpy.dot(
			self._state._moleculeMass[self.containerIndexes],
			self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex]
			)


	def _countsIs(self, values):
		assert (np.size(values) == np.size(self.containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex] = values


	def _countsInc(self, values):
		assert (np.size(values) == np.size(self.containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex] += values


	def _countsDec(self, values):
		assert (np.size(values) == np.size(self.containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex] -= values


class BulkMoleculesView(BulkMoleculesViewBase):
	def __init__(self, *args, **kwargs):
		super(BulkMoleculesView, self).__init__(*args, **kwargs)

		# State references
		self.containerIndexes = self._state.container._namesToIndexes(self._query)


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
		self.containerIndexes = self._state.container._namesToIndexes((self._query,))


	def count(self):
		return self._counts()


	def mass(self):
		return self._mass()


	def countIs(self, value):
		self._countsIs(value)


	def countInc(self, value):
		self._countsInc(value)


	def countDec(self, value):
		self._countsDec(value)

