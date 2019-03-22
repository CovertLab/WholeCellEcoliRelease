"""
Simulation data for external state

This base class includes all data associated with states external to the cells.

	- environment.nutrient_data: a dictionary including the following keys and their values:
		- externalExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- importExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- importConstrainedExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- importUnconstrainedExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- secretionExchangeMolecules: a list of exchange molecules

	- environment.nutrients_time_series: a dictionary of all time series.

	- environment.nutrients_time_series_label: a string specifying the time series
		used for the current simulation.

	Functions:
	----------
	- _buildEnvironment: initializes the environment using conditions and time series
		from raw_data.

	- _getNutrientData: pulls nutrient data from raw_data, saves it in environment.nutrient_data

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.environment import Environment

from reconstruction.ecoli.dataclasses.state import stateFunctions as sf

import re
import numpy as np

class ExternalState(object):
	""" External State """

	def __init__(self, raw_data, sim_data):

		self.environment = Environment(raw_data, sim_data)

		self._buildEnvironment(raw_data, sim_data)

	def _buildEnvironment(self, raw_data, sim_data):

		self.environment.nutrient_data = self._getNutrientData(raw_data)
		self.environment.nutrients_time_series_label = "000000_basal"

		self.environment.nutrients_time_series = {}
		for label in dir(raw_data.condition.timeseries):
			if label.startswith("__"):
				continue

			self.environment.nutrients_time_series[label] = []
			timeseries = getattr(raw_data.condition.timeseries, label)
			for row in timeseries:
				self.environment.nutrients_time_series[label].append((
					row["time"].asNumber(units.s),
					row["nutrients"].encode("utf-8")
					))

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
