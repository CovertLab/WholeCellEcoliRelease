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

		self.uniqueMoleculeDefinitions = None
		self.submass_diff_names = None
		self._submass_diff_name_to_index = None
		self._molecule_ids = None
		self._molecule_masses = None
		self._no_process_assigned_value = None

		self.division_modes = {}


	def initialize(self, sim, sim_data):
		super(UniqueMolecules, self).initialize(sim, sim_data)

		# Used to store information for cell division
		self.uniqueMoleculeDefinitions = deepcopy(sim_data.internal_state.uniqueMolecules.uniqueMoleculeDefinitions)

		# Add the submass difference attributes for processes to operate
		defaultMassAttributes = {}
		self.submass_diff_names = []
		self._submass_diff_name_to_index = {}

		for i, submassName in enumerate(sim_data.submassNameToIndex.viewkeys()):
			massDiffPropertyName = "massDiff_" + submassName
			defaultMassAttributes[massDiffPropertyName] = np.float64
			self.submass_diff_names.append(massDiffPropertyName)
			self._submass_diff_name_to_index[massDiffPropertyName] = i

		for molDef in self.uniqueMoleculeDefinitions.viewvalues():
			molDef.update(defaultMassAttributes)

		self.container = UniqueObjectsContainer(self.uniqueMoleculeDefinitions, self.submass_diff_names)

		# Get ordered list of molecule ids and masses
		self._molecule_ids = tuple(sorted(self.uniqueMoleculeDefinitions.keys()))

		molecule_id_to_mass = {}
		for (id, mass) in izip(
			sim_data.internal_state.uniqueMolecules.uniqueMoleculeMasses["id"],
			sim_data.internal_state.uniqueMolecules.uniqueMoleculeMasses["mass"]
			):
			molecule_id_to_mass[id] = (mass/sim_data.constants.nAvogadro).asNumber(units.fg)

		self._molecule_masses = np.array(
			[molecule_id_to_mass[x] for x in self._molecule_ids]
			)

		self._molecule_id_to_index = {
			x: self._molecule_ids.index(x) for x in self._molecule_ids
			}

		self._no_process_assigned_value = self._nProcesses

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
		processes, all masses are marked as "no_process_assigned".
		"""
		masses = np.zeros(self._masses.shape[1:], np.float64)

		for moleculeId, moleculeMasses in izip(
				self._molecule_ids, self._molecule_masses):
			# Get all molecules of a particular type
			molecules = self.container.objectsInCollection(moleculeId)
			n_molecules = len(molecules)

			if n_molecules == 0:
				continue

			# Add basal masses of the molecule to last row
			masses[self._no_process_assigned_value, :] += moleculeMasses * n_molecules

			# Add additional submasses of the molecule to last row
			massDiffs = molecules.attrsAsStructArray(*self.submass_diff_names).view(
				(np.float64, len(self.submass_diff_names))
				)
			masses[self._no_process_assigned_value, :] += massDiffs.sum(axis=0)

		self._masses[self._preEvolveStateMassIndex, ...] = masses


	def merge(self):
		"""
		Apply all requested changes to the container.
		"""
		self.container.merge()


	def calculatePostEvolveStateMass(self):
		"""
		Computes the summed masses of all unique molecules, after
		evolveState(). If certain process added or deleted molecules, or
		changed the massdiff attributes of molecules, the corresponding change
		of mass is recorded in the self._masses array, with the rows
		corresponding to each process and columns corresponding to each
		submass. This is done to check if each process conserves mass at each
		timestep.
		"""
		# Get mass calculated before evolveState()
		masses = self._masses[self._preEvolveStateMassIndex, ...].copy()

		# Get list of requests
		requests = self.container.get_requests()

		# Loop through all requests
		for req in requests:
			# Apply mass added/subtracted from submass requests
			if req["type"] == "submass":
				process_index = req["source_process_index"]

				for attribute, values in req["added_masses"].viewitems():
					submass_index = self._submass_diff_name_to_index[attribute]
					masses[process_index, submass_index] += values.sum()

			# Delete mass of removed objects
			if req["type"] == "delete":
				process_index = req["source_process_index"]

				# Subtract masses of the deleted molecules themselves
				deleted_masses = self._molecule_masses[req["collection_indexes"], :].sum(axis=0)
				masses[process_index, :] -= deleted_masses

				# Subtract submasses of the deleted molecules
				masses[process_index, :] -= req["deleted_submasses"]

			# Add mass from new objects
			if req["type"] == "new_molecule":
				process_index = req["source_process_index"]

				# Add masses of the added molecules themselves
				collection_index = self._molecule_id_to_index[req["collectionName"]]
				masses_per_molecule = self._molecule_masses[collection_index, :]
				masses[process_index, :] += masses_per_molecule * req["nObjects"]

				# Add submass differences that the molecules were initialized with
				for attribute, values in req["attributes"].viewitems():
					if attribute in self.submass_diff_names:
						submass_index = self._submass_diff_name_to_index[attribute]
						masses[process_index, submass_index] += values.sum()

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

		# self._query must be the name of a unique molecule
		assert isinstance(self._query, basestring)


	def _updateQuery(self):
		# TODO: generalize this logic (both here and in the state)
		# Note: this defaults to a read-only view to the objects
		self._queryResult = self._state.container.objectsInCollection(
			self._query,
			process_index=self._processIndex,
			access=Access.READ_ONLY
			)

		self._totalIs(len(self._queryResult))


	def molecules_read_only(self):
		"""
		Returns a UniqueObjectSet corresponding to the given query with
		read_only access to attributes. The process cannot edit the attributes
		of the molecules or delete molecules with this view.
		"""
		return self._queryResult


	def molecules_read_and_edit(self):
		"""
		Returns a UniqueObjectSet corresponding to the given query with
		read and edit access to attributes. The process cannot delete molecules
		with this view.
		"""
		return self._state.container.objectsInCollection(
			self._query,
			process_index=self._processIndex,
			access=Access.READ_EDIT
			)


	# TODO (ggsun): deprecated alias, should be deleted
	allMolecules = molecules_read_only


	def molecules(self):
		"""
		Returns a UniqueObjectSet corresponding to the given query with full
		read, edit, and delete access.
		"""
		return self._state.container.objectsInCollection(
			self._query,
			process_index=self._processIndex,
			access=Access.READ_EDIT_DELETE
			)


	def moleculeNew(self, **attributes):
		"""
		Adds a single object of the same type as the queried molecule to the
		container with the given attributes.
		"""
		self._state.container.add_new_molecule_request(
			self._query,
			1,
			process_index=self._processIndex,
			**attributes
			)


	def moleculesNew(self, nMolecules, **attributes):
		"""
		Adds nMolecules objects of the same type as the queried molecule to the
		container with the given attributes.
		"""
		self._state.container.add_new_molecule_request(
			self._query,
			nMolecules,
			process_index=self._processIndex,
			**attributes
			)
