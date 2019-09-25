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

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import

import numpy as np

import wholecell.states.external_state
import wholecell.views.view

from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

from wholecell.utils.make_media import Media

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L

ASSERT_POSITIVE_CONCENTRATIONS = True

class NegativeConcentrationError(Exception):
	pass


class LocalEnvironment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):
		self.container = None
		self._moleculeIDs = None
		self._concentrations = None

		self._env_delta_counts = None

		super(LocalEnvironment, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data, timeline):
		super(LocalEnvironment, self).initialize(sim, sim_data, timeline)

		self._processIDs = sim.processes.keys()

		# load target transport reactions from compartment
		boundary_reactions = sim._boundary_reactions
		self.transport_fluxes = {reaction: 0.0 for reaction in boundary_reactions}

		# load constants
		self._nAvogadro = sim_data.constants.nAvogadro

		# make media object
		make_media = Media()

		# if current_timeline_id is specified by a variant in sim_data, look it up in saved_timelines.
		# else, construct the timeline given to initialize
		if sim_data.external_state.environment.current_timeline_id:
			self.current_timeline = sim_data.external_state.environment.saved_timelines[
				sim_data.external_state.environment.current_timeline_id]
		else:
			self.current_timeline = make_media.make_timeline(timeline)

		self.saved_media = sim_data.external_state.environment.saved_media
		self.current_media_id = self.current_timeline[0][1]
		self.current_media = self.saved_media[self.current_media_id]
		self._times = [t[0] for t in self.current_timeline]

		# initialize molecule IDs and concentrations based on initial environment
		self._moleculeIDs = [molecule_id for molecule_id, concentration in self.current_media.iteritems()]
		self._concentrations = np.array([concentration for molecule_id, concentration in self.current_media.iteritems()])
		self._env_delta_counts = dict((molecule_id, 0) for molecule_id in self._moleculeIDs)

		# create bulk container for molecule concentrations. This uses concentrations instead of counts.
		self.container = BulkObjectsContainer(self._moleculeIDs, dtype=np.float64)
		self.container.countsIs(self._concentrations)

		# set the maximum length for a media_id saved to the listener, this is used for padding
		self._media_id_max_length = 25


	def update(self):
		'''update self.current_media_id based on self.current_timeline and self.time'''

		current_index = [i for i, t in enumerate(self._times) if self.time()>=t][-1]

		if self.current_media_id != self.current_timeline[current_index][1]:
			self.current_media_id = self.current_timeline[current_index][1]
			self._concentrations = np.array([concentration for id, concentration in self.saved_media[self.current_media_id].iteritems()])
			self.container.countsIs(self._concentrations)
			print('update media: {}'.format(self.current_media_id))

		if ASSERT_POSITIVE_CONCENTRATIONS and (self._concentrations < 0).any():
			raise NegativeConcentrationError(
					"Negative environment concentration(s) in self._concentrations:\n"
					+ "\n".join("{}".format(self._moleculeIDs[molIndex])
					for molIndex in np.where(self._concentrations < 0)[0]))

	## Functions for multi-scaling interface
	def set_local_environment(self, update):
		# apply environment's concentrations
		concentrations = update['concentrations']
		self._env_delta_counts = dict.fromkeys(self._env_delta_counts, 0)
		for idx, molecule_id in enumerate(self._moleculeIDs):
			self._concentrations[idx] = concentrations[molecule_id]

		# media_id passed from external overwrites the default timeline
		# TODO (eran) -- fix this so that the current_timeline does not need to be re-written
		self.current_media_id = update['media_id']
		self.current_timeline = [(0.0, self.current_media_id)]

		# get transport fluxes passed in from the environment
		self.transport_fluxes = update.get('transport_fluxes', {})

	def get_environment_change(self):
		return self._env_delta_counts

	def accumulate_deltas(self, molecule_ids, counts):
		for molecule_id, count in zip(molecule_ids, counts):
			self._env_delta_counts[molecule_id] += count

	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		tableWriter.writeAttributes(
			# nutrientTimeSeriesLabel = self.current_timeline_id,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			media_id = self.current_media_id.ljust(self._media_id_max_length),
			media_concentrations = self._concentrations,
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
		self._state.accumulate_deltas(molecule_ids, counts)
		return
