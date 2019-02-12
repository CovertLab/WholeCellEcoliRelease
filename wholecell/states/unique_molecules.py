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
import collections
from copy import deepcopy

import numpy as np

import wholecell.states.internal_state
import wholecell.views.view
from wholecell.containers.unique_objects_container import (
	UniqueObjectsContainer, _partition, Access
	)
from wholecell.utils import units


DEFAULT_ATTRIBUTES = {
	"_partitionedProcess":np.int64
	}


class UniqueMolecules(wholecell.states.internal_state.InternalState):
	"""
	UniqueMolecules

	State that tracks unique instances of molecules in the simulation, which
	can have special dynamic attributes.
	"""

	_name = "UniqueMolecules"

	def __init__(self):
		super(UniqueMolecules, self).__init__()

		self.container = None

		self._submassNameToProperty = collections.OrderedDict()

		self.uniqueMoleculeDefinitions = None
		self.submassNameToIndex = None
		self._moleculeIds = None
		self._moleculeMasses = None
		self._unassignedPartitionedValue = None
		self.submass_diff_names = None
		self._mass_changes = None

		self.division_modes = {}


	def initialize(self, sim, sim_data):
		super(UniqueMolecules, self).initialize(sim, sim_data)

		# Used to store information for cell division
		# Should not contain DEFAULT_ATTRIBUTES
		self.uniqueMoleculeDefinitions = deepcopy(sim_data.internal_state.uniqueMolecules.uniqueMoleculeDefinitions)

		# Used to send information to the container
		# Should contain DEFAULT_ATTRIBUTES
		molDefs = sim_data.internal_state.uniqueMolecules.uniqueMoleculeDefinitions.copy()

		defaultAttributes = DEFAULT_ATTRIBUTES.copy()

		self.submassNameToIndex = sim_data.submassNameToIndex

		# Add the submass difference attributes for processes to operate
		defaultMassAttributes = {}
		self.submass_diff_names = []

		for submassName in self.submassNameToIndex.viewkeys():
			massDiffPropertyName = "massDiff_" + submassName
			defaultMassAttributes[massDiffPropertyName] = np.float64
			self._submassNameToProperty[submassName] = massDiffPropertyName
			self.submass_diff_names.append(massDiffPropertyName)

		for molDef in self.uniqueMoleculeDefinitions.viewvalues():
			molDef.update(defaultMassAttributes)

		for molDef in molDefs.viewvalues():
			molDef.update(defaultAttributes)
			molDef.update(defaultMassAttributes)

		self.container = UniqueObjectsContainer(molDefs, self.submass_diff_names)

		self._moleculeIds = sim_data.internal_state.uniqueMolecules.uniqueMoleculeMasses["id"]
		self._moleculeMasses = (
			sim_data.internal_state.uniqueMolecules.uniqueMoleculeMasses["mass"] / sim_data.constants.nAvogadro
			).asNumber(units.fg)[np.argsort(self._moleculeIds)]
		self._moleculeIds = np.sort(self._moleculeIds)

		self._unassignedPartitionedValue = self._nProcesses

		self.division_modes['active_ribosome'] = sim_data.moleculeGroups.unique_molecules_active_ribosome_division
		self.division_modes['domain_index'] = sim_data.moleculeGroups.unique_molecules_domain_index_division
		self.division_modes['binomial'] = sim_data.moleculeGroups.unique_molecules_binomial_division


	def partition(self):
		# Remove any prior partition assignments
		objects = self.container.objects(access=Access.READ_EDIT)
		if len(objects) > 0:
			objects.attrIs(
				_partitionedProcess = self._unassignedPartitionedValue
				)

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
				np.where(partitionedMolecules[:, view._processIndex])[0],
				access=Access.READ_EDIT
				)

			if len(molecules):
				molecules.attrIs(
					_partitionedProcess = view._processIndex,
					)


	def calculatePreEvolveStateMass(self):
		"""
		Computes the summed masses of all unique molecules, prior to
		evolveState(). Since no unique molecules are partitioned to specific
		processes, all masses are marked as "unassigned".
		"""

		if self.simulationStep() == 0:
			# Set everything to the "unassigned" value
			# TODO: consider allowing a default value option for unique objects
			objects = self.container.objects(access=Access.READ_EDIT)

			if len(objects) > 0:
				objects.attrIs(
					_partitionedProcess = self._unassignedPartitionedValue,
					)

		masses = np.zeros(self._masses.shape[1:], np.float64)

		for moleculeId, moleculeMasses in izip(
				self._moleculeIds, self._moleculeMasses):
			# Get all molecules of a particular type
			molecules = self.container.objectsInCollection(moleculeId)
			n_molecules = len(molecules)

			if n_molecules == 0:
				continue

			# Add basal masses of the molecule to last row
			masses[self._unassignedPartitionedValue, :] += moleculeMasses*n_molecules

			# Add additional masses of the molecule to last row
			massDiffs = molecules.attrsAsStructArray(
				*self.submass_diff_names).view(
				(np.float64, len(self.submass_diff_names))
				)
			masses[self._unassignedPartitionedValue, :] += massDiffs.sum(axis=0)

		self._masses[self._preEvolveStateMassIndex, ...] = masses


	def merge(self):
		"""
		Apply all edits made to named attributes of unique objects in the
		container.
		WARNING: This does not check whether conflicting requests were made on
		the same attribute of the same unique molecule.
		"""
		self.container.merge()


	def calculatePostEvolveStateMass(self):
		"""
		Computes the summed masses of all unique molecules, after
		evolveState(). If certain process added or deleted molecules, or
		changed the mass differences of molecules, the corresponding change
		of mass is assigned to the index of the process.
		"""
		# Get mass calculated before evolveState()
		masses = self._masses[self._preEvolveStateMassIndex, ...].copy()

		# Get list of requests
		_, submass_requests, delete_requests, new_molecule_requests = self.container.get_requests()

		for req in submass_requests:
			for attribute, values in req["added_masses"].viewitems():
				submass_index = self.submass_diff_names.index(attribute)
				masses[req["source_process_index"], submass_index] += values.sum()

		for req in delete_requests:
			# Subtract masses of the deleted molecules themselves
			deleted_masses = self._moleculeMasses[req["collection_indexes"], :].sum(axis=0)
			masses[req["source_process_index"], :] -= deleted_masses

			# Subtract deleted submasses
			masses[req["source_process_index"], :] -= req["deleted_submasses"]

		for req in new_molecule_requests:
			# Add masses of the added molecules themselves
			collection_index = np.where(self._moleculeIds == req["collectionName"])[0][0]
			masses_per_molecule = self._moleculeMasses[collection_index, :]
			masses[req["source_process_index"], :] += masses_per_molecule * req["nObjects"]

			# Add submass differences
			for attribute, values in req["attr"].viewitems():
				if attribute in self.submass_diff_names:
					submass_index = self.submass_diff_names.index(attribute)
					masses[req["source_process_index"], submass_index] += values.sum()

		self._masses[self._postEvolveStateMassIndex, ...] = masses

		# Reset request lists in container
		self.container.reset_requests()


	def _calculateMass(self):
		masses = np.zeros(self._masses.shape[1:], np.float64)

		submassDiffNames = self._submassNameToProperty.values() # TODO: cache

		for moleculeId, moleculeMasses in izip(self._moleculeIds, self._moleculeMasses):
			molecules = self.container.objectsInCollection(moleculeId)

			if len(molecules) == 0:
				continue

			processIndexes = molecules.attr("_partitionedProcess")

			countPerProcess = np.bincount(processIndexes, minlength = self._nProcesses + 1)

			masses += np.outer(countPerProcess, moleculeMasses)

			massDiffs = molecules.attrsAsStructArray(*submassDiffNames).view((np.float64, len(submassDiffNames)))

			for processIndex in np.arange(self._nProcesses + 1):
				masses[processIndex, :] += massDiffs[processIndex == processIndexes, :].sum(axis = 0)

		return masses


	def loadSnapshot(self, container):
		"""Copy contents from `container`, which must have the same specifications."""
		self.container.loadSnapshot(container)


	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)


	#def tableAppend(self, tableWriter):
		#self.container.tableAppend(tableWriter)


class UniqueMoleculesView(wholecell.views.view.View):
	_stateID = "UniqueMolecules"

	def __init__(self, *args, **kwargs):
		super(UniqueMoleculesView, self).__init__(*args, **kwargs)

		self._queryResult = None # TODO: store query results with the state

		if isinstance(self._query[0], basestring):
			self._query = list(self._query)
			self._query[0] = [self._query[0]]


	def _updateQuery(self):
		# TODO: generalize this logic (both here and in the state)
		# Note: this defaults to a read-only view to the objects
		self._queryResult = self._state.container.objectsInCollections(
			self._query[0],
			**self._query[1]
			)

		self._totalIs(len(self._queryResult))


	def molecules_read_only(self):
		return self._queryResult


	def molecules_read_and_edit(self):
		return self._state.container.objectsInCollections(
			self._query[0],
			process_index=self._processIndex,
			access=Access.READ_EDIT,
			**self._query[1]
			)


	# TODO (ggsun): deprecated alias, should be deleted
	allMolecules = molecules_read_only


	def molecules(self):
		return self._state.container.objectsInCollections(
			self._query[0],
			process_index=self._processIndex,
			access=Access.READ_EDIT_DELETE,
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
		self._state.container.objectNew(
			moleculeName,
			process_index=self._processIndex,
			_partitionedProcess = self._processIndex,
			**attributes
			)


	def moleculesNew(self, moleculeName, nMolecules, **attributes):
		self._state.container.objectsNew(
			moleculeName,
			nMolecules,
			process_index=self._processIndex,
			_partitionedProcess = self._processIndex,
			**attributes
			)
