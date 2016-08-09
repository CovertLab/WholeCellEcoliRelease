"""
SimulationData for Ecoli

Raw data processed into forms convienent for whole-cell modeling

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""
from __future__ import division

import numpy as np
import collections

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

# Data classes
from reconstruction.ecoli.dataclasses.getterFunctions import getterFunctions
from reconstruction.ecoli.dataclasses.moleculeGroups import moleculeGroups
from reconstruction.ecoli.dataclasses.constants import Constants
from reconstruction.ecoli.dataclasses.state.state import State
from reconstruction.ecoli.dataclasses.process.process import Process
from reconstruction.ecoli.dataclasses.growthRateDependentParameters import Mass, GrowthRateParameters
from reconstruction.ecoli.dataclasses.relation import Relation

from wholecell.utils import units

VERBOSE = False

class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):

		# Doubling time (used in fitting)
		self.doubling_time = None

	def initialize(self, raw_data, basal_expression_condition = "M9 Glucose minus AAs"):

		self._addConditionData(raw_data)
		self.nutrientData = self._getNutrientData(raw_data)
		self.condition = "basal"
		self.nutrientsTimeSeriesLabel = "000000_basal"
		self.doubling_time = self.conditionToDoublingTime[self.condition]

		# TODO: Check that media condition is valid
		self.basal_expression_condition = basal_expression_condition
		# self.envDict, self.externalExchangeMolecules, self.nutrientExchangeMolecules, self.secretionExchangeMolecules = self._addEnvironments(raw_data)

		self._addHardCodedAttributes()

		# Helper functions (have no dependencies)
		self.getter = getterFunctions(raw_data, self)
		self.moleculeGroups = moleculeGroups(raw_data, self)
		self.constants = Constants(raw_data, self)

		# Growth rate dependent parameters are set first
		self.growthRateParameters = GrowthRateParameters(raw_data, self)
		self.mass = Mass(raw_data, self)

		# Data classes (can depend on helper functions)
		# Data classes cannot depend on each other
		self.process = Process(raw_data, self)
		self.state = State(raw_data, self)

		# Relations between data classes (can depend on data classes)
		# Relations cannot depend on each other
		self.relation = Relation(raw_data, self)

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

	def _addEnvironments(self, raw_data):
		externalExchangeMolecules = {}
		nutrientExchangeMolecules = {}
		secretionExchangeMolecules = set()
		envDict = {}
		notEnvList = ["condition_doubling_time", "tf_condition", "condition_defs"]
		environments = [(x, getattr(raw_data.condition, x)) for x in dir(raw_data.condition) if not x.startswith("__") and x not in notEnvList]
		for envName, env in environments:
			externalExchangeMolecules[envName] = set()
			nutrientExchangeMolecules[envName] = set()
			envDict[envName] = collections.deque()
			setpoints = [(float(x.split("_")[-1]), getattr(env, x)) for x in dir(env) if not x.startswith("__")]
			for time, nutrientBounds in setpoints:
				constrainedExchangeMolecules = {}
				unconstrainedExchangeMolecules = []
				for nutrient in nutrientBounds:
					if not np.isnan(nutrient["lower bound"].asNumber()) and not np.isnan(nutrient["upper bound"].asNumber()):
						continue
					elif not np.isnan(nutrient["upper bound"].asNumber()):
						constrainedExchangeMolecules[nutrient["molecule id"]] = nutrient["upper bound"]
						externalExchangeMolecules[envName].add(nutrient["molecule id"])
						nutrientExchangeMolecules[envName].add(nutrient["molecule id"])
					else:
						unconstrainedExchangeMolecules.append(nutrient["molecule id"])
						externalExchangeMolecules[envName].add(nutrient["molecule id"])
						nutrientExchangeMolecules[envName].add(nutrient["molecule id"])

				for secretion in raw_data.secretions:
					if secretion["lower bound"] and secretion["upper bound"]:
						# "non-growth associated maintenance", not included in our metabolic model
						continue

					else:
						externalExchangeMolecules[envName].add(secretion["molecule id"])
						secretionExchangeMolecules.add(secretion["molecule id"])

				D = {
					"constrainedExchangeMolecules": constrainedExchangeMolecules,
					"unconstrainedExchangeMolecules": unconstrainedExchangeMolecules,
					}
				envDict[envName].append((time, D))
			externalExchangeMolecules[envName] = sorted(externalExchangeMolecules[envName])
			nutrientExchangeMolecules[envName] = sorted(nutrientExchangeMolecules[envName])
		secretionExchangeMolecules = sorted(secretionExchangeMolecules)


		return envDict, externalExchangeMolecules, nutrientExchangeMolecules, secretionExchangeMolecules

	def _addConditionData(self, raw_data):
		self.conditionToDoublingTime = dict([(x["condition"].encode("utf-8"), x["doubling time"]) for x in raw_data.condition.condition_defs])
		self.conditionActiveTfs = dict([(x["condition"].encode("utf-8"), x["active TFs"]) for x in raw_data.condition.condition_defs])

		abbrToActiveId = dict([(x["TF"].encode("utf-8"), x["activeId"].encode("utf-8").split(", ")) for x in raw_data.tfIds if len(x["activeId"]) > 0])

		geneIdToRnaId = dict([(x["id"].encode("utf-8"), x["rnaId"].encode("utf-8")) for x in raw_data.genes])

		abbrToRnaId = dict(
			[(x["symbol"].encode("utf-8"), x["rnaId"].encode("utf-8")) for x in raw_data.genes] +
			[(x["name"].encode("utf-8"), geneIdToRnaId[x["geneId"].encode("utf-8")]) for x in raw_data.translationEfficiency if x["geneId"] != "#N/A"]
			)

		self.tfToFC = {}
		notFound = []
		for row in raw_data.foldChanges:
			tf = abbrToActiveId[row["TF"].encode("utf-8")][0]
			try:
				target = abbrToRnaId[row["Target"].encode("utf-8")]
			except KeyError:
				notFound.append(row["Target"].encode("utf-8"))
				continue
			if tf not in self.tfToFC:
				self.tfToFC[tf] = {}
			FC = row["F_avg"]
			if row["Regulation_direct"] < 0:
				FC *= -1.
			FC = 2**FC
			self.tfToFC[tf][target] = FC

		if VERBOSE:
			print "The following target genes listed in foldChanges.tsv have no corresponding entry in genes.tsv:"
			for item in notFound:
				print item


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
				print "Warning: overwriting TF fold change conditions for %s" % tf

			self.tfToActiveInactiveConds[tf]["active genotype perturbations"] = activeGenotype
			self.tfToActiveInactiveConds[tf]["active nutrients"] = activeNutrients
			self.tfToActiveInactiveConds[tf]["inactive genotype perturbations"] = inactiveGenotype
			self.tfToActiveInactiveConds[tf]["inactive nutrients"] = inactiveNutrients

		self.conditions = {}
		for row in raw_data.condition.condition_defs:
			condition = row["condition"].encode("utf-8")
			self.conditions[condition] = {}
			self.conditions[condition]["nutrients"] = row["nutrients"].encode("utf-8")
			self.conditions[condition]["perturbations"] = row["genotype perturbations"]

		for tf in sorted(self.tfToActiveInactiveConds):
			activeCondition = tf + "__active"
			inactiveCondition = tf + "__inactive"
			self.conditions[activeCondition] = {}
			self.conditions[inactiveCondition] = {}
			self.conditions[activeCondition]["nutrients"] = self.tfToActiveInactiveConds[tf]["active nutrients"]
			self.conditions[inactiveCondition]["nutrients"] = self.tfToActiveInactiveConds[tf]["inactive nutrients"]
			self.conditions[activeCondition]["perturbations"] = self.tfToActiveInactiveConds[tf]["active genotype perturbations"]
			self.conditions[inactiveCondition]["perturbations"] = self.tfToActiveInactiveConds[tf]["inactive genotype perturbations"]

		self.nutrientsTimeSeries = {}
		for label in dir(raw_data.condition.timeseries):
			if label.startswith("__"):
				continue

			self.nutrientsTimeSeries[label] = collections.deque()
			timeseries = getattr(raw_data.condition.timeseries, label)
			for row in timeseries:
				self.nutrientsTimeSeries[label].append((
					row["time"].asNumber(units.s),
					row["nutrients"].encode("utf-8")
					))

		self.nutrientToDoublingTime = {}
		for condition in self.conditionToDoublingTime:
			if len(self.conditions[condition]["perturbations"]) > 0:
				continue
			nutrientLabel = self.conditions[condition]["nutrients"]
			if nutrientLabel in self.nutrientToDoublingTime and self.conditionToDoublingTime[condition] != self.nutrientToDoublingTime[nutrientLabel]:
				raise Exception, "Multiple doubling times correspond to the same media conditions"
			self.nutrientToDoublingTime[nutrientLabel] = self.conditionToDoublingTime[condition]


	def _getNutrientData(self, raw_data):

		externalExchangeMolecules = {}
		importExchangeMolecules = {}
		secretionExchangeMolecules = set()
		importConstrainedExchangeMolecules = {}
		importUnconstrainedExchangeMolecules = {}
		nutrientsList = [(x, getattr(raw_data.condition.nutrient, x)) for x in dir(raw_data.condition.nutrient) if not x.startswith("__")]
		for nutrientsName, nutrients in nutrientsList:
			externalExchangeMolecules[nutrientsName] = set()
			importExchangeMolecules[nutrientsName] = set()
			importConstrainedExchangeMolecules[nutrientsName] = {}
			importUnconstrainedExchangeMolecules[nutrientsName] = []
			for nutrient in nutrients:
				if not np.isnan(nutrient["lower bound"].asNumber()) and not np.isnan(nutrient["upper bound"].asNumber()):
					continue
				elif not np.isnan(nutrient["upper bound"].asNumber()):
					importConstrainedExchangeMolecules[nutrientsName][nutrient["molecule id"]] = nutrient["upper bound"]
					externalExchangeMolecules[nutrientsName].add(nutrient["molecule id"])
					importExchangeMolecules[nutrientsName].add(nutrient["molecule id"])
				else:
					importUnconstrainedExchangeMolecules[nutrientsName].append(nutrient["molecule id"])
					externalExchangeMolecules[nutrientsName].add(nutrient["molecule id"])
					importExchangeMolecules[nutrientsName].add(nutrient["molecule id"])

			for secretion in raw_data.secretions:
				if secretion["lower bound"] and secretion["upper bound"]:
					# "non-growth associated maintenance", not included in our metabolic model
					continue

				else:
					externalExchangeMolecules[nutrientsName].add(secretion["molecule id"])
					secretionExchangeMolecules.add(secretion["molecule id"])

			externalExchangeMolecules[nutrientsName] = sorted(externalExchangeMolecules[nutrientsName])
			importExchangeMolecules[nutrientsName] = sorted(importExchangeMolecules[nutrientsName])
		secretionExchangeMolecules = sorted(secretionExchangeMolecules)

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}
