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

		for label in dir(raw_data.condition.timelines):
			if label.startswith("__"):
				continue
			self.environment.nutrients_time_series[label] = []
			timelines = getattr(raw_data.condition.timelines, label)
			for row in timelines:
				self.environment.nutrients_time_series[label].append((
					row["time"].asNumber(units.s),
					row["media"].encode("utf-8"),
					))

		# create a dictionary with all saved environments, including molecules of concentration == 0
		self.environment.environment_dict = {}
		for label in dir(raw_data.condition.media):
			if label.startswith("__"):
				continue
			self.environment.environment_dict[label] = {}

			#initiate all molecules with 0 concentrations
			for row in raw_data.condition.environment_molecules:
				self.environment.environment_dict[label].update({row["molecule id"]: 0})

			# update non-zero concentrations, remove units
			molecule_concentrations = getattr(raw_data.condition.media, label)
			for row in molecule_concentrations:
				self.environment.environment_dict[label].update({row["molecule id"]: row["concentration"].asNumber()})

		# make mapping from external molecule to exchange molecule
		self.environment.env_to_exchange_map = {
			mol["molecule id"]: mol["molecule id"] + mol["exchange molecule location"]
			for mol_index, mol in enumerate(raw_data.condition.environment_molecules)
			}
		self.environment.exchange_to_env_map = {v: k for k, v in self.environment.env_to_exchange_map.viewitems()}

		# make dict with exchange molecules for all saved environments, using env_to_exchange_map
		self.environment.exchange_dict = {}
		for media, concentrations in self.environment.environment_dict.iteritems():
			self.environment.exchange_dict[media] = {
				self.environment.env_to_exchange_map[mol]: conc
				for mol, conc in concentrations.iteritems()
				}

		# initial state based on default nutrient time series
		self.environment.nutrients = self.environment.environment_dict[
			self.environment.nutrients_time_series[self.environment.nutrients_time_series_label][0][1]]
