"""
SimulationData for Ecoli

Raw data processed into forms convenient for whole-cell modeling

"""

from __future__ import absolute_import, division, print_function

import collections

import numpy as np

# Data classes
from reconstruction.ecoli.dataclasses.getter_functions import GetterFunctions
from reconstruction.ecoli.dataclasses.molecule_groups import MoleculeGroups
from reconstruction.ecoli.dataclasses.molecule_ids import MoleculeIds
from reconstruction.ecoli.dataclasses.constants import Constants
from reconstruction.ecoli.dataclasses.common_names import CommonNames
from reconstruction.ecoli.dataclasses.state.internal_state import InternalState
from reconstruction.ecoli.dataclasses.state.external_state import ExternalState
from reconstruction.ecoli.dataclasses.process.process import Process
from reconstruction.ecoli.dataclasses.growth_rate_dependent_parameters import Mass, GrowthRateParameters
from reconstruction.ecoli.dataclasses.relation import Relation
from reconstruction.ecoli.dataclasses.adjustments import Adjustments


VERBOSE = False


class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):

		# Doubling time (used in fitting)
		self.doubling_time = None

	def initialize(self, raw_data, basal_expression_condition="M9 Glucose minus AAs"):
		self._add_condition_data(raw_data)
		self.condition = "basal"
		self.doubling_time = self.condition_to_doubling_time[self.condition]

		# TODO: Check that media condition is valid
		self.basal_expression_condition = basal_expression_condition

		self._add_molecular_weight_keys(raw_data)
		self._add_compartment_keys(raw_data)
		self._add_base_codes(raw_data)

		# General helper functions (have no dependencies)
		self.common_names = CommonNames(raw_data)
		self.constants = Constants(raw_data)
		self.adjustments = Adjustments(raw_data)

		# Reference helper functions (can depend on hard-coded attributes)
		self.molecule_groups = MoleculeGroups(raw_data, self)
		self.molecule_ids = MoleculeIds(raw_data, self)

		# Getter functions (can depend on helper functions and reference classes)
		self.getter = GetterFunctions(raw_data, self)

		# Growth rate dependent parameters are set first
		self.growth_rate_parameters = GrowthRateParameters(raw_data, self)
		self.mass = Mass(raw_data, self)

		# Data classes (can depend on helper and getter functions)
		# Data classes cannot depend on each other
		self.external_state = ExternalState(raw_data, self)
		self.process = Process(raw_data, self)
		self.internal_state = InternalState(raw_data, self)

		# Relations between data classes (can depend on data classes)
		# Relations cannot depend on each other
		self.relation = Relation(raw_data, self)

		self.translation_supply_rate = {}


	def _add_molecular_weight_keys(self, raw_data):
		self.submass_name_to_index = {
			mw_key["submass_name"]: mw_key["index"]
			for mw_key in raw_data.molecular_weight_keys
			}


	def _add_compartment_keys(self, raw_data):
		self.compartment_abbrev_to_index = {
			compartment["abbrev"]: i
			for i, compartment in enumerate(raw_data.compartments)
		}
		self.compartment_id_to_index = {
			compartment["id"]: i
			for i,compartment in enumerate(raw_data.compartments)
		}


	def _add_base_codes(self, raw_data):
		self.amino_acid_code_to_id_ordered = collections.OrderedDict(
			tuple((row["code"], row["id"])
				  for row in raw_data.base_codes.amino_acids))

		self.ntp_code_to_id_ordered = collections.OrderedDict(
			tuple((row["code"], row["id"])
				  for row in raw_data.base_codes.ntp))

		self.dntp_code_to_id_ordered = collections.OrderedDict(
			tuple((row["code"], row["id"])
				  for row in raw_data.base_codes.dntp))


	def _add_condition_data(self, raw_data):
		abbrToActiveId = {x["TF"]: x["activeId"].split(", ") for x in raw_data.transcription_factors if len(x["activeId"]) > 0}
		geneIdToRnaId = {x["id"]: x['rna_id'] for x in raw_data.genes}
		abbrToRnaId = {x["symbol"]: x['rna_id'] for x in raw_data.genes}
		abbrToRnaId.update({
			x["name"]: geneIdToRnaId[x["geneId"]]
			for x in raw_data.translation_efficiency
			if x["geneId"] != "#N/A"})

		self.tf_to_fold_change = {}
		self.tf_to_direction = {}

		removed_fcs = {(row['TF'], row['Target']) for row in raw_data.fold_changes_removed}
		for fc_file in ['fold_changes', 'fold_changes_nca']:
			gene_not_found = set()
			tf_not_found = set()
			for row in getattr(raw_data, fc_file):
				FC = row['log2 FC mean']

				# Skip fold changes that have been removed
				if (row['TF'], row['Target']) in removed_fcs:
					continue

				# Skip fold changes that do not agree with curation
				if row['Regulation_direct'] != '' and row['Regulation_direct'] > 2:
					continue

				# Skip positive autoregulation
				if row['TF'] == row['Target'] and FC > 0:
					continue

				try:
					tf = abbrToActiveId[row['TF']][0]
				except KeyError:
					tf_not_found.add(row['TF'])
					continue

				try:
					target = abbrToRnaId[row['Target']]
				except KeyError:
					gene_not_found.add(row['Target'])
					continue

				if tf not in self.tf_to_fold_change:
					self.tf_to_fold_change[tf] = {}
					self.tf_to_direction[tf] = {}

				self.tf_to_direction[tf][target] = np.sign(FC)
				self.tf_to_fold_change[tf][target] = 2**FC

			if VERBOSE:
				if gene_not_found:
					print(f'The following target genes listed in {fc_file}.tsv'
						' have no corresponding entry in genes.tsv:')
					for item in gene_not_found:
						print(item)

				if tf_not_found:
					print('The following transcription factors listed in'
						f' {fc_file}.tsv have no corresponding active entry in'
						' transcription_factors.tsv:')
					for tf in tf_not_found:
						print(tf)

		self.tf_to_active_inactive_conditions = {}
		for row in raw_data.condition.tf_condition:
			tf = row["active TF"]
			activeGenotype = row["active genotype perturbations"]
			activeNutrients = row["active nutrients"]
			inactiveGenotype = row["inactive genotype perturbations"]
			inactiveNutrients = row["inactive nutrients"]

			if tf not in self.tf_to_active_inactive_conditions:
				self.tf_to_active_inactive_conditions[tf] = {}
			else:
				print("Warning: overwriting TF fold change conditions for %s" % tf)

			self.tf_to_active_inactive_conditions[tf]["active genotype perturbations"] = activeGenotype
			self.tf_to_active_inactive_conditions[tf]["active nutrients"] = activeNutrients
			self.tf_to_active_inactive_conditions[tf]["inactive genotype perturbations"] = inactiveGenotype
			self.tf_to_active_inactive_conditions[tf]["inactive nutrients"] = inactiveNutrients

		# Populate combined conditions data from condition_defs
		self.conditions = {}
		self.condition_to_doubling_time = {}
		self.condition_active_tfs = {}
		self.condition_inactive_tfs = {}
		self.ordered_conditions = []  # order for variant to run
		for row in raw_data.condition.condition_defs:
			condition = row["condition"]
			self.ordered_conditions.append(condition)
			self.conditions[condition] = {}
			self.conditions[condition]["nutrients"] = row["nutrients"]
			self.conditions[condition]["perturbations"] = row["genotype perturbations"]
			self.condition_to_doubling_time[condition] = row['doubling time']
			self.condition_active_tfs[condition] = row['active TFs']
			self.condition_inactive_tfs[condition] = row['inactive TFs']

		# Populate nutrientToDoubling for each set of combined conditions
		self.nutrient_to_doubling_time = {}
		for condition in self.condition_to_doubling_time:
			if len(self.conditions[condition]["perturbations"]) > 0:
				continue
			nutrientLabel = self.conditions[condition]["nutrients"]
			if nutrientLabel in self.nutrient_to_doubling_time and self.condition_to_doubling_time[condition] != self.nutrient_to_doubling_time[nutrientLabel]:
				raise Exception(
					"Multiple doubling times correspond to the same media conditions")
			self.nutrient_to_doubling_time[nutrientLabel] = self.condition_to_doubling_time[condition]

		# Populate conditions and conditionToDboulingTime for active and inactive TF conditions
		basal_dt = self.condition_to_doubling_time['basal']
		for tf in sorted(self.tf_to_active_inactive_conditions):
			for status in ['active', 'inactive']:
				condition = '{}__{}'.format(tf, status)
				nutrients = self.tf_to_active_inactive_conditions[tf]['{} nutrients'.format(status)]
				self.conditions[condition] = {}
				self.conditions[condition]['nutrients'] = nutrients
				self.conditions[condition]['perturbations'] = self.tf_to_active_inactive_conditions[tf]['{} genotype perturbations'.format(status)]
				self.condition_to_doubling_time[condition] = self.nutrient_to_doubling_time.get(nutrients, basal_dt)
