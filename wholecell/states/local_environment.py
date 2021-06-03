#!/usr/bin/env python

"""
External state that represents environmental molecules and conditions.

	Variables:
		- current_media (dict): with {molecule_id: concentration}
		- current_timeline (list): a list of events as tuples with (time, media_id)
		- current_media_id (str): the current media's id
		- current_timeline_id (str): the current timeline's id
		- saved_media (dict): all saved media, keys are media_ids
		- saved_timelines (dict): all saved timelines, keys are timeline_ids
		- _times: a list of all times at which the media shifts

	Functions:
		- update(): updates current_media_id according to saved_timelines

"""

from __future__ import absolute_import, division, print_function

from typing import Dict, Set, Tuple

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
import wholecell.states.external_state
from wholecell.utils import units
import six
from six.moves import zip


COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L

ASSERT_POSITIVE_CONCENTRATIONS = True


class NegativeConcentrationError(Exception):
	pass


class LocalEnvironment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):
		super(LocalEnvironment, self).__init__(*args, **kwargs)

		self.container = None
		self._moleculeIDs = None
		self._env_delta_counts = None

	def initialize(self, sim, sim_data, timeline):
		super(LocalEnvironment, self).initialize(sim, sim_data, timeline)

		# load target transport reactions from compartment
		boundary_reactions = sim._boundary_reactions
		self.transport_fluxes = {reaction: 0.0 for reaction in boundary_reactions}

		# make media object
		make_media = sim_data.external_state.make_media

		# if current_timeline_id is specified by a variant in sim_data, look it up in saved_timelines.
		# else, construct the timeline given to initialize
		if sim_data.external_state.current_timeline_id:
			self.current_timeline = sim_data.external_state.saved_timelines[
				sim_data.external_state.current_timeline_id]
		else:
			self.current_timeline = make_media.make_timeline(timeline)

		self.saved_media = sim_data.external_state.saved_media
		self.current_media_id = self.current_timeline[0][1]
		current_media = self.saved_media[self.current_media_id]
		self._times = [t[0] for t in self.current_timeline]

		# initialize molecule IDs and concentrations based on initial environment
		self._moleculeIDs = [molecule_id for molecule_id, concentration in six.viewitems(current_media)]
		concentrations = np.array([current_media[molecule_id] for molecule_id in self._moleculeIDs])
		self._env_delta_counts = dict((molecule_id, 0) for molecule_id in self._moleculeIDs)

		# create bulk container for molecule concentrations. This uses concentrations instead of counts.
		self.container = BulkObjectsContainer(self._moleculeIDs, dtype=np.float64)
		self.container.countsIs(concentrations)

		# set the maximum length for a media_id saved to the listener, this is used for padding
		self._media_id_max_length = 25

		# Setup for exchange to environment
		self._exchange_data_from_concentrations = sim_data.external_state.exchange_data_from_concentrations
		self.exchange_to_env_map = sim_data.external_state.exchange_to_env_map
		self.import_constraint_threshold = sim_data.external_state.import_constraint_threshold

	def update(self):
		'''update self.current_media_id based on self.current_timeline and self.time'''

		current_index = [i for i, t in enumerate(self._times) if self.time()>=t][-1]

		if self.current_media_id != self.current_timeline[current_index][1]:
			self.current_media_id = self.current_timeline[current_index][1]
			current_media = self.saved_media[self.current_media_id]
			concentrations = np.array([current_media[molecule_id] for molecule_id in self._moleculeIDs])
			self.container.countsIs(concentrations)
			print('update media: {}'.format(self.current_media_id))

		if ASSERT_POSITIVE_CONCENTRATIONS and (self.container.counts() < 0).any():
			raise NegativeConcentrationError(
				"Negative environment concentration(s):\n"
				+ "\n".join("{}".format(self._moleculeIDs[molIndex])
				for molIndex in np.where(self.container.counts() < 0)[0]))

	## Functions for multi-scaling interface
	def set_local_environment(self, update):
		# apply environment's concentrations
		self._env_delta_counts = dict.fromkeys(self._env_delta_counts, 0)
		concentrations = np.array([update['concentrations'][molecule_id] for molecule_id in self._moleculeIDs])
		self.container.countsIs(concentrations)

		# media_id passed from external overwrites the default timeline
		# TODO (eran) -- fix this so that the current_timeline does not need to be re-written
		self.current_media_id = update['media_id']
		self.current_timeline = [(0.0, self.current_media_id)]

		# get transport fluxes passed in from the environment
		self.transport_fluxes = update.get('transport_fluxes', {})

	def get_environment_change(self):
		return self._env_delta_counts

	def accumulate_deltas(self, molecule_ids, counts):
		# TODO: use _env_delta_counts to update?
		for molecule_id, count in zip(molecule_ids, counts):
			self._env_delta_counts[molecule_id] += count

	def get_exchange_data(self):
		# type: () -> Tuple[Set[str], Dict[str, units.Unum]]
		current_concentrations = dict(zip(self._moleculeIDs, self.container.counts()))
		exchange_data = self._exchange_data_from_concentrations(current_concentrations)
		unconstrained = exchange_data['importUnconstrainedExchangeMolecules']
		constrained = exchange_data['importConstrainedExchangeMolecules']
		return unconstrained, constrained

	def get_import_molecules(self):
		unconstrained, constrained = self.get_exchange_data()
		imports = set(unconstrained) | set(constrained)
		return imports

	def molecule_exchange(self, exchange_molecules, counts):
		'''
		Convert exchange molecules to environmental molecules using mapping
		and updates deltas in the environment.

		Args:
			exchange_molecules (tuple[str]): internal molecule ID that exchanges
				with the environment
			counts (np.ndarray[int]): change in counts for all exchange molecules
		'''

		molecule_ids = [self.exchange_to_env_map[m] for m in exchange_molecules]
		self.accumulate_deltas(molecule_ids, counts)

	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		tableWriter.writeAttributes(
			# nutrientTimeSeriesLabel = self.current_timeline_id,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			media_id = self.current_media_id.ljust(self._media_id_max_length),
			media_concentrations = self.container.counts(),
			)


class EnvironmentView(object):
	_stateID = 'Environment'

	def __init__(self, state, process, query): # weight, priority, coupling id, option to not evaluate the query
		self._state = state
		self._state.viewAdd(self)
		self._processId = process.name()
		self._processIndex = process._processIndex
		self._query = query

		assert len(set(self._query)) == len(self._query), "Environment views cannot contain duplicate entries"
		self._containerIndexes = self._state.container._namesToIndexes(self._query)

	def _dataSize(self):
		return len(self._query)

	def totalConcentrations(self):
		return self._state.container.counts()[self._containerIndexes]

	def import_present(self):
		return self.totalConcentrations() > self._state.import_constraint_threshold
