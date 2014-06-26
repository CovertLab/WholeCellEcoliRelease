"""
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueObjectsContainer object, which 
creates and manages the structured arrays in memory.
"""

from __future__ import division

from itertools import izip

import numpy as np
import tables

import wholecell.states.state
import wholecell.views.view
from wholecell.containers.unique_objects_container import UniqueObjectsContainer, _partition

DEFAULT_ATTRIBUTES = {
	"_partitionedProcess":np.int64
	}

class UniqueMolecules(wholecell.states.state.State):
	"""
	UniqueMolecules

	State that tracks unique instances of molecules in the simulation, which 
	can have special dynamic attributes.
	"""

	_name = "UniqueMolecules"

	def __init__(self, *args, **kwargs):
		self.container = None

		self._submassNameToProperty = {}

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		molDefs = kb.uniqueMoleculeDefinitions.copy()

		defaultAttributes = DEFAULT_ATTRIBUTES.copy()

		self.submassNameToIndex = kb.submassNameToIndex

		# Add the submass difference attributes for processes to operate
		for submassName in self.submassNameToIndex.viewkeys():
			massDiffPropertyName = "massDiff_" + submassName
			defaultAttributes[massDiffPropertyName] = np.float64
			self._submassNameToProperty[submassName] = massDiffPropertyName

		for molDef in molDefs.viewvalues():
			molDef.update(defaultAttributes)

		self.container = UniqueObjectsContainer(molDefs)

		self._moleculeIds = kb.uniqueMoleculeMasses["moleculeId"]
		self._moleculeMasses = kb.uniqueMoleculeMasses["mass"].to("fg/mole").magnitude / kb.nAvogadro.magnitude

		self._unassignedPartitionedValue = self._nProcesses


	def partition(self):
		# Set the correct time for saving purposes
		self.container.timeIs(self.timeStep())

		# Remove any prior partition assignments
		self.container.objects().attrIs(_partitionedProcess = self._unassignedPartitionedValue)
		
		# Gather requests
		nMolecules = self.container._globalReference.size
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

		# Don't calculate on non-requesting views
		doCalculatePartition = (requestNumberVector > 0)

		objectRequestsArray[:, ~doCalculatePartition] = False

		partitionedMolecules = np.zeros((nMolecules, self._nProcesses), np.bool)

		# Grant non-overlapping requests all of the relevant molecules
		overlappingRequests = np.dot(objectRequestsArray.T, objectRequestsArray)

		overlappingRequests[np.identity(nViews, np.bool)] = False

		# TODO: ignore overlapping requests with a process

		noOverlap = ~overlappingRequests.any(0)

		for viewIndex in np.where(noOverlap)[0]:
			# NOTE: there may be a way to vectorize this
			partitionedMolecules[:,self._views[viewIndex]._processIndex] |= objectRequestsArray[:, viewIndex]

		doCalculatePartition[noOverlap] = False

		if doCalculatePartition.any():
			partitionedMolecules |= _partition(
				objectRequestsArray[:, doCalculatePartition],
				requestNumberVector[doCalculatePartition],
				requestProcessArray[doCalculatePartition, :],
				self.randomState
				)

		for view in self._views:
			molecules = self.container.objectsByGlobalIndex(
				np.where(partitionedMolecules[:, view._processIndex])[0]
				)

			if len(molecules):
				molecules.attrIs(_partitionedProcess = view._processIndex)

		# Compute partitioned masses

		masses = np.zeros(self._masses.shape[1:], np.float64)

		submassDiffNames = self._submassNameToProperty.values() # TODO: cache

		for moleculeId, moleculeMasses in izip(self._moleculeIds, self._moleculeMasses):
			molecules = self.container.objectsInCollection(moleculeId)

			processIndexes = molecules.attr("_partitionedProcess")

			countPerProcess = np.bincount(processIndexes, minlength = self._nProcesses + 1)

			masses += np.outer(countPerProcess, moleculeMasses)

			massDiffs = molecules.attrsAsStructArray(*submassDiffNames).view((np.float64, len(submassDiffNames)))

			masses[processIndexes, :] += massDiffs

		self._masses[self._preEvolveStateMassIndex, ...] = masses


	def merge(self):
		# Operations are performed directly on the container, so there is no
		# "merge" operation needed

		# Compute partitioned masses

		masses = np.zeros(self._masses.shape[1:], np.float64)

		submassDiffNames = self._submassNameToProperty.values() # TODO: cache

		for moleculeId, moleculeMasses in izip(self._moleculeIds, self._moleculeMasses):
			molecules = self.container.objectsInCollection(moleculeId)

			processIndexes = molecules.attr("_partitionedProcess")

			countPerProcess = np.bincount(processIndexes, minlength = self._nProcesses + 1)

			masses += np.outer(countPerProcess, moleculeMasses)

			massDiffs = molecules.attrsAsStructArray(*submassDiffNames).view((np.float64, len(submassDiffNames)))

			masses[processIndexes, :] += massDiffs

		self._masses[self._postEvolveStateMassIndex, ...] = masses


	# TODO: refactor mass calculations as a whole
	def mass(self):
		totalMass = 0

		submassDiffNames = self._submassNameToProperty.values()
		
		for moleculeId, masses in izip(self._moleculeIds, self._moleculeMasses):
			molecules = self.container.objectsInCollection(moleculeId)

			nMolecules = len(molecules)

			if nMolecules == 0:
				continue

			totalMass += (
				masses.sum() * nMolecules
				+ molecules.attrsAsStructArray(
					*submassDiffNames
					).astype(np.float64).sum()
				)

		return totalMass


	def massByType(self, typeKey):
		totalMass = 0

		submassIndex = self.submassNameToIndex[typeKey]
		submassDiffName = self._submassNameToProperty[typeKey]
		
		for moleculeId, masses in izip(self._moleculeIds, self._moleculeMasses):
			molecules = self.container.objectsInCollection(moleculeId)

			nMolecules = len(molecules)

			if nMolecules == 0:
				continue

			totalMass += (
				masses[submassIndex] * nMolecules
				+ molecules.attr(
					submassDiffName
					).sum()
				)

		return totalMass


	def massByCompartment(self, compartment):
		return 0


	def pytablesCreate(self, h5file, expectedRows):
		# self.container.pytablesCreate(h5file)
		pass


	def pytablesAppend(self, h5file):
		# self.container.pytablesAppend(h5file)
		pass


	def pytablesLoad(self, h5file, timePoint):
		# self.container.pytablesLoad(h5file, timePoint)
		raise Exception("Unique molecules saving disabled for now")


class UniqueMoleculesView(wholecell.views.view.View):
	__slots__ = ("_state", "_processId", "_processIndex", "_query", "_totalCount", "_requestedCount", "_queryResult")
	_stateID = "UniqueMolecules"

	def __init__(self, *args, **kwargs):
		super(UniqueMoleculesView, self).__init__(*args, **kwargs)

		self._queryResult = None # TODO: store query results with the state


	def _updateQuery(self):
		# TODO: generalize this logic (both here and in the state)

		self._queryResult = self._state.container.objectsInCollection(
			self._query[0],
			**self._query[1]
			)

		self._totalIs(len(self._queryResult))

	def allMolecules(self):
		return self._queryResult

	def molecules(self):
		return self._state.container.objectsInCollection(
			self._query[0],
			_partitionedProcess = ("==", self._processIndex),
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
		return self._state.container.objectNew(
			moleculeName,
			_partitionedProcess = self._processIndex,
			**attributes
			)


	def moleculesNew(self, moleculeName, nMolecules, **attributes):
		return self._state.container.objectsNew(
			moleculeName,
			nMolecules,
			_partitionedProcess = self._processIndex,
			**attributes
			)
