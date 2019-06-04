"""
Simulation data for external state

This base class includes all data associated with states external to the cells.
Initializes the environment using conditions and time series from raw_data.

	- environment.saved_timelines: a dictionary of all timelines.
	- environment.current_timeline_id: a string specifying the timelines
		used for the current simulation.
	- environment.current_media: a dictionary of molecules (keys) and
		their concentrations (values).
	- environment.saved_media: a dictionary of all media, each entry
		itself a dictionary molecules (keys) and their concentrations (values).

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

from wholecell.utils import units
from reconstruction.ecoli.dataclasses.state.environment import Environment
from wholecell.utils.make_media import Media


class ExternalState(object):
	""" External State """

	def __init__(self, raw_data, sim_data):
		self.environment = Environment(raw_data, sim_data)

		# make media object
		make_media = Media()

		# create a dictionary with all saved timelines
		self.environment.saved_timelines = {}
		for row in raw_data.condition.timelines_def:
			timeline_id = row["timeline"]
			timeline_str = row["events"]
			new_timeline = make_media.make_timeline(timeline_str)
			self.environment.saved_timelines[timeline_id] = new_timeline

		# set default current_timeline_id to None, this can be overwritten by the timelines variant
		self.environment.current_timeline_id = None

		# make a dictionary with all media conditions specified by media_recipes
		make_media = Media()
		self.environment.saved_media = make_media.make_saved_media()

		# make mapping from external molecule to exchange molecule
		self.environment.env_to_exchange_map = {
			mol["molecule id"]: mol["molecule id"] + mol["exchange molecule location"]
			for mol_index, mol in enumerate(raw_data.condition.environment_molecules)
			}
		self.environment.exchange_to_env_map = {v: k for k, v in self.environment.env_to_exchange_map.viewitems()}

		# make dict with exchange molecules for all saved environments, using env_to_exchange_map
		self.environment.exchange_dict = {}
		for media, concentrations in self.environment.saved_media.iteritems():
			self.environment.exchange_dict[media] = {
				self.environment.env_to_exchange_map[mol]: conc
				for mol, conc in concentrations.iteritems()
				}
