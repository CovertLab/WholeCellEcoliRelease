"""
Simulation data for external state

This base class includes all data associated with states external to the cells.
Initializes the environment using conditions and time series from raw_data.

	- saved_timelines: a dictionary of all timelines.
	- current_timeline_id: a string specifying the timelines
		used for the current simulation.
	- current_media: a dictionary of molecules (keys) and
		their concentrations (values).
	- saved_media: a dictionary of all media, each entry
		itself a dictionary molecules (keys) and their concentrations (values).

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

from typing import Any, Dict

import numpy as np

from wholecell.utils import units
from wholecell.utils.make_media import Media
import six


# threshold (units.mmol / units.L) separates concentrations that are import constrained with
# max flux = 0 from unconstrained molecules.
IMPORT_CONSTRAINT_THRESHOLD =  1e-5


class ExternalState(object):
	""" External State """

	def __init__(self, raw_data, sim_data):
		# make media object
		self.make_media = Media(raw_data)

		self.carbon_sources = sim_data.moleculeGroups.carbon_sources
		self._initialize_environment(raw_data)
		self.all_external_exchange_molecules = self._get_all_external_exchange_molecules(raw_data)
		self.secretion_exchange_molecules = self._get_secretion_exchange_molecules(raw_data)


	def _get_all_external_exchange_molecules(self, raw_data):
		'''
		Returns:
			list[str]: all external exchange molecules
		'''
		externalExchangeData = []
		# initiate all molecules with 0 concentrations
		for row in raw_data.condition.environment_molecules:
			externalExchangeData.append(row["molecule id"] + row["exchange molecule location"])

		return externalExchangeData

	def _get_secretion_exchange_molecules(self, raw_data):
		'''
		Returns:
			set[str]: all secretion exchange molecules
		'''
		secretionExchangeMolecules = []
		for secretion in raw_data.secretions:
			if secretion["lower bound"] and secretion["upper bound"]:
				# "non-growth associated maintenance", not included in our metabolic model
				continue
			else:
				secretionExchangeMolecules.append(secretion["molecule id"])

		return set(secretionExchangeMolecules)

	def _initialize_environment(self, raw_data):
		self.import_constraint_threshold = IMPORT_CONSTRAINT_THRESHOLD

		# create a dictionary with all saved timelines
		self.saved_timelines = {}
		for row in raw_data.condition.timelines_def:
			timeline_id = row["timeline"]
			timeline_str = row["events"]
			new_timeline = self.make_media.make_timeline(timeline_str)
			self.saved_timelines[timeline_id] = new_timeline

		# set default current_timeline_id to None, this can be overwritten by the timelines variant
		self.current_timeline_id = None

		# make a dictionary with all media conditions specified by media_recipes
		self.saved_media = self.make_media.make_saved_media()

		# make mapping from external molecule to exchange molecule
		self.env_to_exchange_map = {
			mol["molecule id"]: mol["molecule id"] + mol["exchange molecule location"]
			for mol_index, mol in enumerate(raw_data.condition.environment_molecules)
			}
		self.exchange_to_env_map = {v: k for k, v in six.viewitems(self.env_to_exchange_map)}

		# make dict with exchange molecules for all saved environments, using env_to_exchange_map
		self.exchange_dict = {}
		for media, concentrations in six.viewitems(self.saved_media):
			self.exchange_dict[media] = {
				self.env_to_exchange_map[mol]: conc
				for mol, conc in six.viewitems(concentrations)
				}

	def exchange_data_from_concentrations(self, molecules):
		# type: (Dict[str, float]) -> Dict[str, Any]
		'''
		Update importExchangeMolecules for FBA based on current nutrient concentrations.
		This provides a simple type of transport to accommodate changing nutrient
		concentrations in the environment. Transport is modeled as a binary switch:
		When there is a high concentrations of environment nutrients, transporters
		are unconstrained and nutrients are transported as needed by metabolism.
		When concentrations fall below the threshold, that nutrient's transport
		is constrained to max flux of 0.

		Args:
			molecules: external molecules (no location tag) with external concentration,
				concentration can be inf

		Returns dict with the following keys:
			externalExchangeMolecules (set[str]): all exchange molecules (with
				location tag), includes both import and secretion exchanged molecules
			importExchangeMolecules (set[str]): molecules (with location tag) that
				can be imported from the environment into the cell
			importConstrainedExchangeMolecules (dict[str, float with mol/mass/time units]):
				constrained molecules (with location tag) with upper bound flux constraints
			importUnconstrainedExchangeMolecules (set[str]): exchange molecules
				(with location tag) that do not have an upper bound on their flux
			secretionExchangeMolecules (set[str]): molecules (with location tag)
				that can be secreted by the cell into the environment
		'''

		externalExchangeMolecules = set()
		importExchangeMolecules = set()
		secretionExchangeMolecules = self.secretion_exchange_molecules

		oxygen_id = 'OXYGEN-MOLECULE[p]'

		exchange_molecules = {self.env_to_exchange_map[mol]: conc for mol, conc in six.viewitems(molecules)}

		# Unconstrained uptake if greater than import threshold
		importUnconstrainedExchangeMolecules = {molecule_id
			for molecule_id, concentration in exchange_molecules.items()
			if concentration >= self.import_constraint_threshold}
		importExchangeMolecules.update(importUnconstrainedExchangeMolecules)
		externalExchangeMolecules.update(importUnconstrainedExchangeMolecules)

		# TODO: functionalize limits based on concentrations of transporters and environment
		# Limit carbon uptake if present depending on the presence of oxygen
		importConstrainedExchangeMolecules = {}
		for carbon_source_id in self.carbon_sources:
			if carbon_source_id in importUnconstrainedExchangeMolecules:
				if oxygen_id in importUnconstrainedExchangeMolecules:
					importConstrainedExchangeMolecules[carbon_source_id] = 20. * (units.mmol / units.g / units.h)
				else:
					importConstrainedExchangeMolecules[carbon_source_id] = 100. * (units.mmol / units.g / units.h)
				importUnconstrainedExchangeMolecules.remove(carbon_source_id)

		externalExchangeMolecules.update(secretionExchangeMolecules)

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}

	def exchange_data_from_media(self, media_label):
		'''
		Returns:
			dict: exchange_data for a media_label saved in exchange_data_dict.
		'''

		concentrations = self.saved_media[media_label]
		return self.exchange_data_from_concentrations(concentrations)

	def get_import_constraints(self, unconstrained, constrained, units):
		'''
		Returns:
			unconstrained_molecules (list[bool]): the indices of all
				importUnconstrainedExchangeMolecules in
				self.all_external_exchange_molecules are true, the rest as false
			constrained_molecules (list[bool]): the indices of all
				importConstrainedExchangeMolecules in
				self.all_external_exchange_molecules are true, the rest as false
			constraints (list[float]): uptake constraints for each molecule
				that is constrained, nan for no constraint
		'''

		# molecules from all_external_exchange_molecules set to 'true' if they are current importExchangeMolecules.
		unconstrained_molecules = [
			molecule_id in unconstrained
			for molecule_id in self.all_external_exchange_molecules
			]

		# molecules from all_external_exchange_molecules set to 'true' if they are current importConstrainedExchangeMolecules.
		constrained_molecules = [
			molecule_id in constrained
			for molecule_id in self.all_external_exchange_molecules
			]

		constraints = [
			constrained.get(molecule_id, np.nan * units).asNumber(units)
			for molecule_id in self.all_external_exchange_molecules
			]

		return unconstrained_molecules, constrained_molecules, constraints
