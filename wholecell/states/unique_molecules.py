"""
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueObjectsContainer object, which
creates and manages the structured arrays in memory.
"""

from __future__ import absolute_import, division, print_function

from copy import deepcopy

import numpy as np
import six
from six.moves import zip

import wholecell.states.internal_state
import wholecell.views.view
from wholecell.containers.unique_objects_container import UniqueObjectsContainer
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

		self._submass_diff_names = None
		self._submass_diff_name_to_index = None
		self._molecule_ids = None
		self._molecule_masses = None
		self._pre_evolve_state_mass_index = None
		self._requests = []

		self.uniqueMoleculeDefinitions = None
		self.division_mode = {}


	def initialize(self, sim, sim_data):
		super(UniqueMolecules, self).initialize(sim, sim_data)

		# Used to store information for cell division
		self.uniqueMoleculeDefinitions = deepcopy(
			sim_data.internal_state.uniqueMolecules.uniqueMoleculeDefinitions)

		# Add the submass difference attributes for processes to operate
		defaultMassAttributes = {}
		self._submass_diff_names = []
		self._submass_diff_name_to_index = {}

		for submass_name, index in sim_data.submass_name_to_index.items():
			massDiffPropertyName = "massDiff_" + submass_name
			defaultMassAttributes[massDiffPropertyName] = np.float64
			self._submass_diff_names.append(massDiffPropertyName)
			self._submass_diff_name_to_index[massDiffPropertyName] = index

		for molDef in six.viewvalues(self.uniqueMoleculeDefinitions):
			molDef.update(defaultMassAttributes)

		self.container = UniqueObjectsContainer(
			self.uniqueMoleculeDefinitions, submass_diff_names=self._submass_diff_names)

		# Get ordered list of molecule ids and masses
		self._molecule_ids = tuple(sorted(self.uniqueMoleculeDefinitions.keys()))

		molecule_id_to_mass = {}
		uniqueMoleculeMasses = sim_data.internal_state.uniqueMolecules.uniqueMoleculeMasses
		for (id_, mass) in zip(
			uniqueMoleculeMasses["id"], uniqueMoleculeMasses["mass"]
			):
			molecule_id_to_mass[id_] = (mass/sim_data.constants.nAvogadro).asNumber(units.fg)

		self._molecule_masses = np.array(
			[molecule_id_to_mass[x] for x in self._molecule_ids]
			)

		self._molecule_id_to_index = {
			x: self._molecule_ids.index(x) for x in self._molecule_ids
			}

		# Total mass of molecules before evolveState goes into the last row
		# of the masses array
		self._pre_evolve_state_mass_index = self._nProcesses

		self.division_mode['active_ribosome'] = sim_data.moleculeGroups.unique_molecules_active_ribosome_division
		self.division_mode['RNA'] = sim_data.moleculeGroups.unique_molecules_RNA_division
		self.division_mode['domain_index'] = sim_data.moleculeGroups.unique_molecules_domain_index_division
		self.division_mode['chromosomal_segment'] = sim_data.moleculeGroups.unique_molecules_chromosomal_segment_division


	def partition(self, processes):
		"""
		Unique molecules are not partitioned.
		"""
		pass


	def merge(self, processes):
		"""
		Apply all requested changes to the container. The returned list of
		requests is used to compute mass differences for each process.
		"""
		self._requests = self.container.merge()

		# Initialize mass differences array
		process_mass_diffs = np.zeros_like(self._process_mass_diffs)

		# Loop through all requests
		for req in self._requests:
			# Apply mass added/subtracted from submass requests
			if req["type"] == "submass":
				process_index = req["process_index"]

				for attribute, values in six.viewitems(req["added_masses"]):
					submass_index = self._submass_diff_name_to_index[attribute]
					process_mass_diffs[process_index, submass_index] += values.sum()

			# Delete mass of removed objects
			if req["type"] == "delete":
				process_index = req["process_index"]

				# Subtract masses of the deleted molecules themselves
				deleted_masses = self._molecule_masses[req["collection_indexes"], :].sum(axis=0)
				process_mass_diffs[process_index, :] -= deleted_masses

				# Subtract submasses of the deleted molecules
				process_mass_diffs[process_index, :] -= req["deleted_submasses"]

			# Add mass from new objects
			if req["type"] == "new_molecule":
				process_index = req["process_index"]

				# Add masses of the added molecules themselves
				collection_index = self._molecule_id_to_index[req["collectionName"]]
				masses_per_molecule = self._molecule_masses[collection_index, :]
				process_mass_diffs[process_index, :] += masses_per_molecule * req["nObjects"]

				# Add submass differences that the molecules were initialized with
				for attribute, values in six.viewitems(req["attributes"]):
					if attribute in self._submass_diff_names:
						submass_index = self._submass_diff_name_to_index[attribute]
						process_mass_diffs[process_index, submass_index] += values.sum()

		self._process_mass_diffs += process_mass_diffs


	def calculateMass(self):
		"""
		Computes the summed masses of all unique molecules, including both the
		basal mass and the added submasses of each molecule.
		"""
		masses = np.zeros_like(self._masses)

		for moleculeId, moleculeMasses in zip(
				self._molecule_ids, self._molecule_masses):
			# Get all molecules of a particular type
			molecules = self.container.objectsInCollection(moleculeId)
			n_molecules = len(molecules)

			if n_molecules == 0:
				continue

			# Add basal masses of the molecule
			masses += moleculeMasses * n_molecules

			# Add additional submasses of the molecule
			massDiffs = molecules.attrsAsStructArray(*self._submass_diff_names).view(
				(np.float64, len(self._submass_diff_names))
				)
			masses += massDiffs.sum(axis=0)

		self._masses = masses


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
		self.cached_attributes = {}

		# self._query must be the name of a unique molecule
		assert isinstance(self._query, six.string_types)

	def _updateQuery(self):
		# TODO: generalize this logic (both here and in the state)
		# Note: this defaults to a read-only view to the objects
		self._queryResult = self._state.container.objectsInCollection(
			self._query,
			process_index=self._processIndex,
			access=()
			)
		self._flush_cached_attributes()

		self._totalIs(len(self._queryResult))

	def _flush_cached_attributes(self):
		"""
		Removes all of the cached attribute arrays.
		"""
		self.cached_attributes = {}

	def request_access(self, access):
		"""
		Requests access permissions required to edit or delete the attributes
		of the unique molecules being viewed. Argument should be a tuple of
		requested access permissions e.g. (Access.EDIT, Access.DELETE)
		"""
		self._queryResult.set_access_level(access=access)

	# Wrappers for reading or manipulating queried molecules
	def attr(self, attribute):
		if attribute not in self.cached_attributes:
			self.cached_attributes[attribute] = self._queryResult.attr(
				attribute)
		return self.cached_attributes[attribute]

	def attrs(self, *attributes):
		all_attrs = []
		for attribute in attributes:
			if attribute not in self.cached_attributes:
				self.cached_attributes[attribute] = self._queryResult.attr(
					attribute)
			all_attrs.append(self.cached_attributes[attribute])

		return tuple(all_attrs)

	def attrIs(self, **attributes):
		self._flush_cached_attributes()
		self._queryResult.attrIs(**attributes)

	def add_submass_by_name(self, submass_name, delta_mass):
		self._queryResult.add_submass_by_name(submass_name, delta_mass)

	def add_submass_by_array(self, delta_mass):
		self._queryResult.add_submass_by_array(delta_mass)

	def delByIndexes(self, indexes):
		self._flush_cached_attributes()
		self._queryResult.delByIndexes(indexes)

	def moleculesNew(self, nMolecules, **attributes):
		"""
		Adds nMolecules objects of the same type as the queried molecule to the
		container with the given attributes.
		"""
		self._flush_cached_attributes()

		unique_indexes = self._state.container.add_request(
			type="new_molecule",
			collectionName=self._query,
			nObjects=nMolecules,
			process_index=self._processIndex,
			attributes=attributes
			)

		return unique_indexes
