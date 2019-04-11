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
from environment.condition.make_media import Media

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

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

		# create a dictionary with all media conditions specified by media_recipes
		self.environment.saved_media = {}
		for row in raw_data.condition.media_recipes:
			new_media_id = row["media id"]
			base_id = row["base media"]
			added_media_id = row["added media"]
			ingredient_ids = row["ingredients"]
			base_media = make_media.stock_media[base_id]

			if added_media_id:
				added_media = make_media.stock_media[added_media_id]
				base_vol = row["base media volume"]
				added_vol = row["added media volume"]
				new_media = make_media.combine_media(base_media, base_vol, added_media, added_vol)
			elif ingredient_ids:
				base_vol = row.get("base media volume", 0 * units.L)
				added_weight = row.get("ingredients weight", None)
				added_counts = row.get("ingredients counts", None)
				added_vol = row.get("ingredients volume")  # the row is a list with units.L, even an empty list is read.

				ingredients = {ingred_id: {} for ingred_id in ingredient_ids}
				for index, ingred_id in enumerate(ingredient_ids):
					if added_weight:
						ingredients[ingred_id]['weight'] = added_weight[index]
					if added_counts:
						ingredients[ingred_id]['counts'] = added_counts[index]
					if added_vol:
						ingredients[ingred_id]['volume'] = added_vol[index]
					else:
						ingredients[ingred_id]['volume'] = 0 * units.L

				new_media = make_media.add_ingredients(base_media, base_vol, ingredients)
			else:
				new_media = base_media

			# remove concentration units, setting at CONC_UNITS
			unitless_new_media = {mol: conc.asNumber(CONC_UNITS) for mol, conc in new_media.iteritems()}
			self.environment.saved_media[new_media_id] = unitless_new_media

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
