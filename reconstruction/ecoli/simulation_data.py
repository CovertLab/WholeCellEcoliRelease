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
		self._add_hard_coded_attributes()

		# General helper functions (have no dependencies)
		self.common_names = CommonNames(raw_data)
		self.constants = Constants(raw_data)

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


	def _add_hard_coded_attributes(self):
		self.amino_acid_code_to_id_ordered = collections.OrderedDict((
			("A", "L-ALPHA-ALANINE[c]"), ("R", "ARG[c]"), ("N", "ASN[c]"), ("D", "L-ASPARTATE[c]"),
			("C", "CYS[c]"), ("E", "GLT[c]"), ("Q", "GLN[c]"), ("G", "GLY[c]"),
			("H", "HIS[c]"), ("I", "ILE[c]"), ("L", "LEU[c]"), ("K", "LYS[c]"),
			("M", "MET[c]"), ("F", "PHE[c]"), ("P", "PRO[c]"), ("S", "SER[c]"),
			("T", "THR[c]"), ("W", "TRP[c]"), ("Y", "TYR[c]"), ("U", "L-SELENOCYSTEINE[c]"),
			("V", "VAL[c]")
			))

		self.ntp_code_to_id_ordered = collections.OrderedDict((
			("A", "ATP[c]"), ("C", "CTP[c]"), ("G", "GTP[c]"), ("U", "UTP[c]")
			))

		self.dntp_code_to_id_ordered = collections.OrderedDict((
			("A", "DATP[c]"), ("C", "DCTP[c]"), ("G", "DGTP[c]"), ("T", "TTP[c]")
			))


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
		notFound = []
		for row in raw_data.fold_changes:
			# Skip fold changes that do not agree with curation
			if np.abs(row['Regulation_direct']) > 2:
				continue

			tf = abbrToActiveId[row["TF"]][0]
			try:
				target = abbrToRnaId[row["Target"]]
			except KeyError:
				notFound.append(row["Target"])
				continue
			if tf not in self.tf_to_fold_change:
				self.tf_to_fold_change[tf] = {}
				self.tf_to_direction[tf] = {}
			FC = row["F_avg"]
			if row["Regulation_direct"] < 0:
				FC *= -1.
				self.tf_to_direction[tf][target] = -1
			else:
				self.tf_to_direction[tf][target] = 1
			FC = 2**FC
			self.tf_to_fold_change[tf][target] = FC

		if VERBOSE:
			print("The following target genes listed in fold_changes.tsv have no corresponding entry in genes.tsv:")
			for item in notFound:
				print(item)

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
		self.ordered_conditions = []  # order for variant to run
		for row in raw_data.condition.condition_defs:
			condition = row["condition"]
			self.ordered_conditions.append(condition)
			self.conditions[condition] = {}
			self.conditions[condition]["nutrients"] = row["nutrients"]
			self.conditions[condition]["perturbations"] = row["genotype perturbations"]
			self.condition_to_doubling_time[condition] = row['doubling time']
			self.condition_active_tfs[condition] = row['active TFs']

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
