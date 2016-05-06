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
from unum import Unum

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

class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):

		# Doubling time (used in fitting)
		self.doubling_time = None

	def initialize(self, doubling_time, raw_data, expression_condition = "M9 Glucose minus AAs", environment = "wildtype"):

		if type(doubling_time) != Unum:
			raise Exception("Doubling time is not a Unum object!")
		self.doubling_time = doubling_time

		# TODO: Check that media condition is valid
		self.expression_condition = expression_condition
		self.envDict, self.externalExchangeMolecules, self.nutrientExchangeMolecules = self._addEnvironments(raw_data)
		self.environment = environment

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
		envDict = {}
		environments = [(x, getattr(raw_data.environment, x)) for x in dir(raw_data.environment) if not x.startswith("__")]
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

				D = {
					"constrainedExchangeMolecules": constrainedExchangeMolecules,
					"unconstrainedExchangeMolecules": unconstrainedExchangeMolecules,
					}
				envDict[envName].append((time, D))
			externalExchangeMolecules[envName] = sorted(externalExchangeMolecules[envName])
			nutrientExchangeMolecules[envName] = sorted(nutrientExchangeMolecules[envName])


		return envDict, externalExchangeMolecules, nutrientExchangeMolecules
