#!/usr/bin/env python

"""
External state that represents environmental molecules and conditions.
	- nutrients_time_series: a list of tuples that include time and nutrients in
		which shifts occur.
	- nutrients: a string specifying the current nutrient condition.
	- times: a list of all times at which the nutrients shift.
	Functions:
	----------
	- update: updates nutrients according to nutrients_time_series
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np

import wholecell.states.external_state
import wholecell.views.view

from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L

ASSERT_POSITIVE_CONCENTRATIONS = True

class NegativeConcentrationError(Exception):
	pass


class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):
		self.container = None
		self._moleculeIDs = None
		self._concentrations = None

		self._environment_deltas = None

		super(Environment, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self._processIDs = sim.processes.keys()

		# load constants
		self._nAvogadro = sim_data.constants.nAvogadro

		# environment time series data
		self.environment_dict = sim_data.external_state.environment.environment_dict
		self.nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label

		# get current nutrients label
		self.current_time_series = sim_data.external_state.environment.nutrients_time_series[self.nutrients_time_series_label]
		self.nutrients = self.current_time_series[0][1]
		self._times = [t[0] for t in self.current_time_series]

		# initialize molecule IDs and concentrations based on initial environment
		self._moleculeIDs = [molecule_id for molecule_id, concentration in self.environment_dict[self.nutrients].iteritems()]
		self._concentrations = np.array([concentration.asNumber() for molecule_id, concentration in self.environment_dict[self.nutrients].iteritems()])

		# create bulk container for molecule concentrations. This uses concentrations instead of counts.
		self.container = BulkObjectsContainer(self._moleculeIDs, dtype=np.float64)
		self.container.countsIs(self._concentrations)

		# the length of the longest nutrients name, for padding in nutrients listener
		self._nutrients_name_max_length = max([len(t[1]) for t in self.current_time_series])


	def update(self):
		current_index = [i for i, t in enumerate(self._times) if self.time()>=t][-1]

		# update nutrients based on nutrient_time_series. This updates the concentrations,
		# and also the nutrients label is used in polypeptide_elongation to find
		# a ribosomeElongationRate in ribosomeElongationRateDict
		if self.nutrients != self.current_time_series[current_index][1]:
			self.nutrients = self.current_time_series[current_index][1]
			self._concentrations = np.array([concentration.asNumber() for id, concentration in self.environment_dict[self.nutrients].iteritems()])
			self.container.countsIs(self._concentrations)

		if ASSERT_POSITIVE_CONCENTRATIONS and (self._concentrations < 0).any():
			raise NegativeConcentrationError(
					"Negative environment concentration(s) in self._concentrations:\n"
					+ "\n".join(
					"{}".format(
						self._moleculeIDs[molIndex],
						)
					for molIndex in np.where(self._concentrations < 0)[0]
					)
				)

	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		tableWriter.writeAttributes(
			nutrientTimeSeriesLabel = self.nutrients_time_series_label,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			nutrientCondition = self.nutrients.ljust(self._nutrients_name_max_length),
			nutrientConcentrations = self._concentrations,
			)

class EnvironmentViewBase(object):
	_stateID = 'Environment'

	def __init__(self, state, process, query): # weight, priority, coupling id, option to not evaluate the query
		self._state = state
		self._state.viewAdd(self)
		self._processId = process.name()
		self._processIndex = process._processIndex
		self._query = query
		self._concentrations = np.zeros(self._dataSize(), np.float64) # number of objects that satisfy the query


	# Interface to State
	def _updateQuery(self):
		self._totalIs(self._state.container._counts[self._containerIndexes])


	def _totalIs(self, value):
		self._concentrations[:] = value


	def _countsInc(self, counts):
		return


	# Interface to Process
	def _totalConcentrations(self):
		return np.array(self._state._concentrations)[self._containerIndexes].copy()



class EnvironmentView(EnvironmentViewBase):
	def __init__(self, *args, **kwargs):
		super(EnvironmentView, self).__init__(*args, **kwargs)

		# State references
		assert len(set(self._query)) == len(self._query), "Environment views cannot contain duplicate entries"
		self._containerIndexes = self._state.container._namesToIndexes(self._query)


	def _dataSize(self):
		return len(self._query)


	def totalConcentrations(self):
		return self._totalConcentrations()


	def countsInc(self, molecule_ids, counts):
		self._state._environment_deltas = counts
		#TODO (Eran) save deltas for external environment dict(zip(molecule_ids, counts))
		#TODO (Eran) deltas size varies because of changing importExchange.  This will need to be fixed for a listener to save these

		return
