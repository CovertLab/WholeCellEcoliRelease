'''
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueObjectsContainer object, which 
creates and manages the structured arrays in memory.
'''

from __future__ import division

import itertools

import numpy as np
import tables

import wholecell.states.state
import wholecell.views.view
from wholecell.containers.unique_objects_container import UniqueObjectsContainer, _partition


DEFAULT_ATTRIBUTES = {
	'massDiffMetabolite':np.float,
	'massDiffRna':np.float,
	'massDiffProtein':np.float,
	}


class UniqueMolecules(wholecell.states.state.State):
	'''
	UniqueMolecules

	State that tracks unique instances of molecules in the simulation, which 
	can have special dynamic attributes.
	'''

	_name = "UniqueMolecules"

	def __init__(self, *args, **kwargs):
		self.time = None

		self.container = None

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		molDefs = kb.uniqueMoleculeDefinitions.copy()

		for molDef in molDefs.viewvalues():
			molDef.update(DEFAULT_ATTRIBUTES)

		self.container = UniqueObjectsContainer(molDefs)

		self._masses = kb.uniqueMoleculeMasses


	def partition(self):
		# Set the correct time for saving purposes
		self.container.timeIs(self.timeStep())

		# Clear out any deleted entries to make room for new molecules
		self.container.flushDeleted()
		
		# Gather requests
		nMolecules = self.container._collections[self.container._globalRefIndex].size
		nViews = len(self._views)

		objectRequestsArray = np.zeros((nMolecules, nViews), np.bool)
		requestNumberVector = np.zeros(nViews, np.int64)
		requestProcessArray = np.zeros((nViews, self._nProcesses), np.bool)

		for viewIndex, view in enumerate(self._views):
			objectRequestsArray[view._queryResult._globalIndexes, viewIndex] = True

			requestNumberVector[viewIndex] = view._request()

			requestProcessArray[viewIndex, view._processIndex] = True

		# TODO: move this logic to the _partition function
		if requestNumberVector.sum() == 0:
			return

		partitionedMolecules = _partition(objectRequestsArray,
			requestNumberVector, requestProcessArray, self.randStream)

		for view in self._views:
			molecules = self.container.objectsByGlobalIndex(
				np.where(partitionedMolecules[:, view._processIndex])[0]
				)

			for molecule in molecules:
				molecule.attrIs(_partitionedProcess = view._processIndex + 1)
				# "0", being the default, is reserved for unpartitioned molecules


	# TODO: refactor mass calculations as a whole
	def mass(self):
		# TODO: rework this so it's a faster operation (a dot product)

		totalMass = 0
		
		for entry in self._masses:
			moleculeId = entry['moleculeId']
			massMetabolite = entry['massMetabolite']
			massRna = entry['massRna']
			massProtein = entry['massProtein']

			molecules = self.container.objectsInCollection(moleculeId)

			if len(molecules) == 0:
				continue

			for molecule in molecules:
				totalMass += massMetabolite
				totalMass += massRna
				totalMass += massProtein

				totalMass += molecule.attr('massDiffMetabolite')
				totalMass += molecule.attr('massDiffRna')
				totalMass += molecule.attr('massDiffProtein')

		return totalMass


	def massByType(self, typeKey):
		# TODO: rework this so it's a faster operation (a dot product)

		if typeKey in ['rrnas', 'water']:
			return 0

		submassKey = {
			'metabolites':'massMetabolite',
			'rnas':'massRna',
			'proteins':'massProtein',
			}[typeKey]

		submassDiffKey = {
			'metabolites':'massDiffMetabolite',
			'rnas':'massDiffRna',
			'proteins':'massDiffProtein',
			}[typeKey]

		totalMass = 0
		
		for entry in self._masses:
			moleculeId = entry['moleculeId']
			mass = entry[submassKey]

			molecules = self.container.objectsInCollection(moleculeId)

			if len(molecules) == 0:
				continue

			for molecule in molecules:
				totalMass += mass

				totalMass += molecule.attr(submassDiffKey)

		return totalMass


	def massByCompartment(self, compartment):
		return 0


	def pytablesCreate(self, h5file, expectedRows):
		self.container.pytablesCreate(h5file)


	def pytablesAppend(self, h5file):
		self.container.pytablesAppend(h5file)


	def pytablesLoad(self, h5file, timePoint):
		self.container.pytablesLoad(h5file, timePoint)


class UniqueMoleculesView(wholecell.views.view.View):
	_stateID = 'UniqueMolecules'

	def __init__(self, *args, **kwargs):
		super(UniqueMoleculesView, self).__init__(*args, **kwargs)

		self._queryResult = None


	def _updateQuery(self):
		# TODO: generalize this logic (both here and in the state)

		self._queryResult = self._state.container.objectsInCollection(
			self._query[0],
			**self._query[1]
			)

		self._totalIs(len(self._queryResult))


	def molecules(self):
		return self._state.container.objectsInCollection(
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
