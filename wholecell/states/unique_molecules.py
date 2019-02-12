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
	UniqueObjectsContainer, Access
	)
from wholecell.utils import units


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
		self.uniqueMoleculeDefinitions = deepcopy(sim_data.internal_state.uniqueMolecules.uniqueMoleculeDefinitions)

		# Used to send information to the container
		molDefs = sim_data.internal_state.uniqueMolecules.uniqueMoleculeDefinitions.copy()

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
		"""
		Unique molecules are not partitioned.
		"""
		pass


	def calculatePreEvolveStateMass(self):
		"""
		Computes the summed masses of all unique molecules, prior to
		evolveState(). Since no unique molecules are partitioned to specific
		processes, all masses are marked as "unassigned".
		"""
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
			process_index=self._processIndex,
			access=Access.READ_ONLY,
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
			**attributes
			)


	def moleculesNew(self, moleculeName, nMolecules, **attributes):
		self._state.container.objectsNew(
			moleculeName,
			nMolecules,
			process_index=self._processIndex,
			**attributes
			)
