"""
SimulationData for Ecoli

Raw data processed into forms convenient for whole-cell modeling

"""

from __future__ import annotations

import collections

import numpy as np
import scipy

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
from wholecell.utils.fitting import normalize


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
		self.pPromoterBound = {}


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
		gene_id_to_rna_id = {
			gene['id']: gene['rna_ids'][0] for gene in raw_data.genes}
		gene_symbol_to_rna_id = {
			gene['symbol']: gene['rna_ids'][0] for gene in raw_data.genes}
		gene_symbol_to_rna_id.update({
			x["name"]: gene_id_to_rna_id[x["geneId"]]
			for x in raw_data.translation_efficiency
			if x["geneId"] != "#N/A"})

		rna_ids_with_coordinates = {
			gene['rna_ids'][0] for gene in raw_data.genes
			if gene['left_end_pos'] is not None and gene['right_end_pos'] is not None}

		self.tf_to_fold_change = {}
		self.tf_to_direction = {}

		for fc_file in ['fold_changes', 'fold_changes_nca']:
			gene_not_found = set()
			tf_not_found = set()
			gene_location_not_specified = set()

			for row in getattr(raw_data, fc_file):
				FC = row['log2 FC mean']

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
					target = gene_symbol_to_rna_id[row['Target']]
				except KeyError:
					gene_not_found.add(row['Target'])
					continue

				if target not in rna_ids_with_coordinates:
					gene_location_not_specified.add(row['Target'])
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

				if gene_location_not_specified:
					print(f'The following target genes listed in {fc_file}.tsv'
						  ' have no chromosomal location specified in'
						  ' genes.tsv:')
					for item in gene_location_not_specified:
						print(item)

		self.tf_to_active_inactive_conditions = {}
		for row in raw_data.condition.tf_condition:
			tf = row["active TF"]

			if tf not in self.tf_to_fold_change:
				continue

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

	def calculate_ppgpp_expression(self, condition: str):
		"""
		Calculates the expected expression of RNA based on ppGpp regulation
		in a given condition and the expected transcription factor effects in
		that condition.

		Relies on other values that are calculated in the fitting process so
		should only be called after the parca has been run.

		Args:
			condition: label for the desired condition to calculate the average
				expression for (eg. 'basal', 'with_aa', etc)
		"""

		ppgpp = self.growth_rate_parameters.get_ppGpp_conc(
			self.condition_to_doubling_time[condition])
		delta_prob = self.process.transcription_regulation.get_delta_prob_matrix(ppgpp=True)
		p_promoter_bound = np.array([
			self.pPromoterBound[condition][tf]
			for tf in self.process.transcription_regulation.tf_ids
			])
		delta = delta_prob @ p_promoter_bound
		prob, factor = self.process.transcription.synth_prob_from_ppgpp(
			ppgpp, self.process.replication.get_average_copy_number)
		rna_expression = prob * (1 + delta) / factor

		# For cases with no basal ppGpp expression, assume the delta prob is the
		# same as without ppGpp control
		mask = prob == 0
		rna_expression[mask] = delta[mask] / factor[mask]

		rna_expression[rna_expression < 0] = 0
		return normalize(rna_expression)
