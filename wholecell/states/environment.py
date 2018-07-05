#!/usr/bin/env python

"""
External state that represents environmental molecules and conditions.

	- nutrient_data: a dictionary including the following keys and their values:
		- externalExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- importExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- importConstrainedExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- importUnconstrainedExchangeMolecules: a dictionary of all the nutrient condition names, with a list of molecules as their values.
		- secretionExchangeMolecules: a list of exchange molecules

	- nutrients_time_series: a list of tuples that include time and nutrients in which shifts occur.

	- nutrients: a string specifying the current nutrient condition.

	- times: a list of all times at which the nutrients shift.

	Functions:
	----------
	- update: updates nutrients according to nutrients_time_series

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np

import wholecell.states.external_state


class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self.nutrient_data = sim_data.external_state.environment.nutrient_data
		self.nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label
		self.nutrients_time_series = sim_data.external_state.environment.nutrients_time_series[
			self.nutrients_time_series_label
			]
		self.nutrients = self.nutrients_time_series[0][1]
		self.times = [t[0] for t in self.nutrients_time_series]

		# save the length of the longest nutrients name, for padding names in listener
		self.nutrients_name_max_length = len(max([t[1] for t in self.nutrients_time_series], key=len))


	def update(self):
		current_index = [i for i, t in enumerate(self.times) if self.time()>=t][-1]
		self.nutrients = self.nutrients_time_series[current_index][1]


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			nutrientTimeSeriesLabel = self.nutrients_time_series_label,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			nutrients = self.nutrients.ljust(self.nutrients_name_max_length),
			)
