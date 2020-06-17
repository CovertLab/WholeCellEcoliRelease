"""
SimulationData for Ecoli

Raw data processed into forms convenient for whole-cell modeling

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import collections

import numpy as np

# Data classes
from reconstruction.ecoli.dataclasses.getterFunctions import getterFunctions
from reconstruction.ecoli.dataclasses.moleculeGroups import MoleculeGroups
from reconstruction.ecoli.dataclasses.moleculeIds import MoleculeIds
from reconstruction.ecoli.dataclasses.constants import Constants
from reconstruction.ecoli.dataclasses.common_names import CommonNames
from reconstruction.ecoli.dataclasses.state.internal_state import InternalState
from reconstruction.ecoli.dataclasses.state.external_state import ExternalState
from reconstruction.ecoli.dataclasses.process.process import Process
from reconstruction.ecoli.dataclasses.growthRateDependentParameters import Mass, GrowthRateParameters
from reconstruction.ecoli.dataclasses.relation import Relation


VERBOSE = False


class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):

		# Doubling time (used in fitting)
		self.doubling_time = None

	def initialize(self, raw_data, basal_expression_condition = "M9 Glucose minus AAs"):
		self._addConditionData(raw_data)
		self.condition = "basal"
		self.doubling_time = self.conditionToDoublingTime[self.condition]

		# TODO: Check that media condition is valid
		self.basal_expression_condition = basal_expression_condition

		self._addHardCodedAttributes()

		# Helper functions (have no dependencies)
		self.getter = getterFunctions(raw_data, self)
		self.moleculeGroups = MoleculeGroups(raw_data, self)
		self.moleculeIds = MoleculeIds(raw_data, self)
		self.constants = Constants(raw_data, self)

		# Growth rate dependent parameters are set first
		self.growthRateParameters = GrowthRateParameters(raw_data, self)
		self.mass = Mass(raw_data, self)

		# Data classes (can depend on helper functions)
		# Data classes cannot depend on each other
		self.external_state = ExternalState(raw_data, self)
		self.process = Process(raw_data, self)
		self.internal_state = InternalState(raw_data, self)

		# Relations between data classes (can depend on data classes)
		# Relations cannot depend on each other
		self.relation = Relation(raw_data, self)

		self.translationSupplyRate = {}

		# Mappings from element IDs to common names
		self.common_names = CommonNames(raw_data, self)


	def _addHardCodedAttributes(self):
		self.molecular_weight_keys = [
			'23srRNA',
			'16srRNA',
			'5srRNA',
			'tRNA',
			'mRNA',
			'miscRNA',
			'protein',
			'metabolite',
			'water',
			'DNA',
			'RNA' # nonspecific RNA
			]

		self.molecular_weight_order = collections.OrderedDict([
			(key, index) for index, key in enumerate(self.molecular_weight_keys)
			])

		self.submassNameToIndex = self.molecular_weight_order

		self.amino_acid_1_to_3_ordered = collections.OrderedDict((
			("A", "L-ALPHA-ALANINE[c]"), ("R", "ARG[c]"), ("N", "ASN[c]"), ("D", "L-ASPARTATE[c]"),
			("C", "CYS[c]"), ("E", "GLT[c]"), ("Q", "GLN[c]"), ("G", "GLY[c]"),
			("H", "HIS[c]"), ("I", "ILE[c]"), ("L", "LEU[c]"), ("K", "LYS[c]"),
			("M", "MET[c]"), ("F", "PHE[c]"), ("P", "PRO[c]"), ("S", "SER[c]"),
			("T", "THR[c]"), ("W", "TRP[c]"), ("Y", "TYR[c]"), ("U", "L-SELENOCYSTEINE[c]"),
			("V", "VAL[c]")
			))

		self.dNtpOrder = ["A", "C", "G", "T"]


	def _addConditionData(self, raw_data):
		abbrToActiveId = dict([(x["TF"].encode("utf-8"), x["activeId"].encode("utf-8").split(", ")) for x in raw_data.tfIds if len(x["activeId"]) > 0])
		geneIdToRnaId = dict([(x["id"].encode("utf-8"), x["rnaId"].encode("utf-8")) for x in raw_data.genes])
		abbrToRnaId = dict(
			[(x["symbol"].encode("utf-8"), x["rnaId"].encode("utf-8")) for x in raw_data.genes] +
			[(x["name"].encode("utf-8"), geneIdToRnaId[x["geneId"].encode("utf-8")]) for x in raw_data.translationEfficiency if x["geneId"] != "#N/A"]
			)

		self.tfToFC = {}
		self.tfToDirection = {}
		notFound = []
		for row in raw_data.foldChanges:
			# Skip fold changes that do not agree with curation
			if np.abs(row['Regulation_direct']) > 2:
				continue

			tf = abbrToActiveId[row["TF"].encode("utf-8")][0]
			try:
				target = abbrToRnaId[row["Target"].encode("utf-8")]
			except KeyError:
				notFound.append(row["Target"].encode("utf-8"))
				continue
			if tf not in self.tfToFC:
				self.tfToFC[tf] = {}
				self.tfToDirection[tf] = {}
			FC = row["F_avg"]
			if row["Regulation_direct"] < 0:
				FC *= -1.
				self.tfToDirection[tf][target] = -1
			else:
				self.tfToDirection[tf][target] = 1
			FC = 2**FC
			self.tfToFC[tf][target] = FC

		if VERBOSE:
			print("The following target genes listed in foldChanges.tsv have no corresponding entry in genes.tsv:")
			for item in notFound:
				print(item)

		self.tfToActiveInactiveConds = {}
		for row in raw_data.condition.tf_condition:
			tf = row["active TF"].encode("utf-8")
			activeGenotype = row["active genotype perturbations"]
			activeNutrients = row["active nutrients"].encode("utf-8")
			inactiveGenotype = row["inactive genotype perturbations"]
			inactiveNutrients = row["inactive nutrients"].encode("utf-8")

			if tf not in self.tfToActiveInactiveConds:
				self.tfToActiveInactiveConds[tf] = {}
			else:
				print("Warning: overwriting TF fold change conditions for %s" % tf)

			self.tfToActiveInactiveConds[tf]["active genotype perturbations"] = activeGenotype
			self.tfToActiveInactiveConds[tf]["active nutrients"] = activeNutrients
			self.tfToActiveInactiveConds[tf]["inactive genotype perturbations"] = inactiveGenotype
			self.tfToActiveInactiveConds[tf]["inactive nutrients"] = inactiveNutrients

		# Populate combined conditions data from condition_defs
		self.conditions = {}
		self.conditionToDoublingTime = {}
		self.conditionActiveTfs = {}
		self.ordered_conditions = []  # order for variant to run
		for row in raw_data.condition.condition_defs:
			condition = row["condition"].encode("utf-8")
			self.ordered_conditions.append(condition)
			self.conditions[condition] = {}
			self.conditions[condition]["nutrients"] = row["nutrients"].encode("utf-8")
			self.conditions[condition]["perturbations"] = row["genotype perturbations"]
			self.conditionToDoublingTime[condition] = row['doubling time']
			self.conditionActiveTfs[condition] = row['active TFs']

		# Populate nutrientToDoubling for each set of combined conditions
		self.nutrientToDoublingTime = {}
		for condition in self.conditionToDoublingTime:
			if len(self.conditions[condition]["perturbations"]) > 0:
				continue
			nutrientLabel = self.conditions[condition]["nutrients"]
			if nutrientLabel in self.nutrientToDoublingTime and self.conditionToDoublingTime[condition] != self.nutrientToDoublingTime[nutrientLabel]:
				raise Exception(
					"Multiple doubling times correspond to the same media conditions")
			self.nutrientToDoublingTime[nutrientLabel] = self.conditionToDoublingTime[condition]

		# Populate conditions and conditionToDboulingTime for active and inactive TF conditions
		basal_dt = self.conditionToDoublingTime['basal']
		for tf in sorted(self.tfToActiveInactiveConds):
			for status in ['active', 'inactive']:
				condition = '{}__{}'.format(tf, status)
				nutrients = self.tfToActiveInactiveConds[tf]['{} nutrients'.format(status)]
				self.conditions[condition] = {}
				self.conditions[condition]['nutrients'] = nutrients
				self.conditions[condition]['perturbations'] = self.tfToActiveInactiveConds[tf]['{} genotype perturbations'.format(status)]
				self.conditionToDoublingTime[condition] = self.nutrientToDoublingTime.get(nutrients, basal_dt)
