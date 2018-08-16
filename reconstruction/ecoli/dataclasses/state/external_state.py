"""
Simulation data for external state

This base class includes all data associated with states external to the cells.
Initializes the environment using conditions and time series from raw_data.

	- environment.nutrients_time_series: a dictionary of all time series.
	- environment.nutrients_time_series_label: a string specifying the time series
		used for the current simulation.
	- environment.nutrients: a dictionary of environmental nutrients (keys) and
		their concentrations (values).
	- environment.environment_dict: a dictionary of all environments, each one
		itself a dictionary nutrients (keys) and their concentrations (values).

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

		# default parameters
		self.environment.nutrients_time_series_label = "000000_basal"

		# create a dictionary of all nutrient time series
		self.environment.nutrients_time_series = {}
		for label in dir(raw_data.condition.timeseries):
			if label.startswith("__"):
				continue
			self.environment.nutrients_time_series[label] = []
			timeseries = getattr(raw_data.condition.timeseries, label)
			for row in timeseries:
				self.environment.nutrients_time_series[label].append((
					row["time"].asNumber(units.s),
					row["nutrients"].encode("utf-8"),
					))

		# create a dictionary with all saved environments, including molecules of concentration == 0
		self.environment.environment_dict = {}
		for label in dir(raw_data.condition.environment):
			if label.startswith("__"):
				continue
			self.environment.environment_dict[label] = {}

			#initiate all molecules with 0 concentrations
			for row in raw_data.condition.environment_molecules:
				self.environment.environment_dict[label].update({row["molecule id"]: 0 * (units.mmol / units.L)})

			# update non-zero concentrations
			molecule_concentrations = getattr(raw_data.condition.environment, label)
			for row in molecule_concentrations:
				self.environment.environment_dict[label].update({row["molecule id"]: row["concentration"]})

		# initial state based on default nutrient time series
		self.environment.nutrients = self.environment.environment_dict[
			self.environment.nutrients_time_series[self.environment.nutrients_time_series_label][0][1]]
