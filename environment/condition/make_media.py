'''
Functions for making media

'''

from __future__ import absolute_import, division, print_function

from wholecell.utils import units

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

INF = float("inf")
NEG_INF = float("-inf")

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS


class AddIngredientsError(Exception):
	pass

class Media(object):
	'''
	A media object is a factory for making new media by either combining
	two saved media at different volumes (with self.combine_media(), or
	adding ingredients to a saved media (with self.add_ingredients()).
	Ingredients can be added by either specifying their weight (in grams)
	or the counts (in mmol) in addition to the volume. The new media dicts
	are returned to the caller, and are not saved in this object. A media
	object holds dicts about stock media in ```self.stock_media``` and the
	formula weight of environmental molecules in
	```self.environment_molecules_fw```, which is needed for mixing in
	ingredients at weights.
	'''

	def __init__(self):

		raw_data = KnowledgeBaseEcoli()

		# get dicts from knowledge base
		self.environment_molecules_fw = self._get_environment_molecules_fw(raw_data)
		self.stock_media = self._get_stock_media(raw_data)

	def _get_environment_molecules_fw(self, raw_data):
		'''get formula weight (units.g / units.mol) for all environmental molecules'''

		environment_molecules_fw = {}
		for row in raw_data.condition.environment_molecules:
			mol = row["molecule id"]
			fw = row["formula weight"]
			if fw == 'None':
				environment_molecules_fw[mol] = None
			else:
				environment_molecules_fw[mol] = float(fw) * (units.g / units.mol)

		return environment_molecules_fw

	def _get_stock_media(self, raw_data):
		'''load all stock media'''

		stock_media = {}
		for label in vars(raw_data.condition.media):
			# initiate all molecules with 0 concentrations
			stock_media[label] = {
				row["molecule id"]: 0.0 * CONC_UNITS
				for row in raw_data.condition.environment_molecules}

			# get non-zero concentrations (assuming units.mmol / units.L)
			molecule_concentrations = getattr(raw_data.condition.media, label)

			environment_non_zero_dict = {
				row["molecule id"]: row["concentration"]
				for row in molecule_concentrations}

			# update saved_media with non zero concentrations
			stock_media[label].update(environment_non_zero_dict)

		return stock_media

	def combine_media(self, base_media, base_media_volume, mix_media, mix_media_volume):
		'''
		Combines two medias and returns a new media

		Args:
			base_media, mix_media (dict): dicts with {molecule_id: concentration}
			base_media_volume, mix_media_volume (unum): the volumes of base_media and mix_media (floats) with a volume units (i.e. units.L)

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		# intialize new_media
		new_media = {mol_id: 0.0 * CONC_UNITS for mol_id, conc in base_media.iteritems()}

		# get new_media volume
		new_volume = base_media_volume + mix_media_volume

		for mol_id, base_conc in base_media.iteritems():
			mix_conc = mix_media[mol_id]

			if base_conc.asNumber() == INF or mix_conc.asNumber() == INF:
				new_media[mol_id] = INF * CONC_UNITS
			else:
				base_counts = base_conc * base_media_volume
				mix_counts = mix_conc * mix_media_volume
				new_counts = base_counts + mix_counts
				new_conc = new_counts / new_volume

				# update media
				new_media[mol_id] = new_conc

		return new_media

	def add_ingredients(self, base_media, base_media_volume, ingredients):
		'''
		Combines ingredients to existing media to make new media.

		Args:
			base_media (dict): {molecule_id: concentrations}
			base_media_volume:
			ingredients (dict): keys are ingredient ids, values are dicts with weight, counts, volume.
				Only one of weights (in g) or counts (in mmol) is needed; if both are specified, it will use weight.
				If weight or counts is Infinity, the new concentration is set to inf. If the weight or counts is -Infinity,
				the new concentration is set to 0.
				Example format of ingredients:
					{mol_id_1: {'weight': 1.78 * units.g, 'volume': 0.025 * units.L),
					mol_id_2: {'counts': 0.2 * units.mmol, 'volume': 0.1 * units.L),
					}

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		# intialize new_media
		new_media = {mol_id: 0.0 * CONC_UNITS for mol_id, conc in base_media.iteritems()}

		# get new_media volume
		ingredients_volume = 0 * VOLUME_UNITS
		for quantities in ingredients.itervalues():
			ingredients_volume += quantities['volume']
		new_volume = base_media_volume + ingredients_volume

		# get new_media concentrations from mixing ingredients
		for mol_id, base_conc in base_media.iteritems():

			if mol_id in ingredients:
				base_counts = base_conc * base_media_volume
				quantities = ingredients[mol_id]
				weight = quantities.get('weight', None)
				mix_counts = quantities.get('counts', None)

				# calculate mix_counts from weight.
				# This will overwrite added counts if those were specified
				if weight is not None:
					if weight.asNumber() == INF:
						mix_counts = INF * COUNTS_UNITS
					elif weight.asNumber() == NEG_INF:
						mix_counts = NEG_INF * COUNTS_UNITS
					elif weight.asNumber() >= 0:
						if self.environment_molecules_fw[mol_id] is not None:
							fw = self.environment_molecules_fw[mol_id]
							mix_counts = weight / fw
						else:
							raise AddIngredientsError(
								"No fw defined for {} in environment_molecules.tsv".format(mol_id)
							)
					else:
						raise AddIngredientsError(
							"Negative weight given for {}".format(mol_id)
						)
				elif mix_counts is None:
					raise AddIngredientsError(
						"No added added weight or counts for {}".format(mol_id)
					)

				# get new concentration
				# make infinite concentration of ingredient if mix_counts is Infinity
				if mix_counts.asNumber() == INF:
					new_media[mol_id] = INF * CONC_UNITS

				# remove ingredient from media if mix_counts is -Infinity
				# this will override infinite concentrations in base_media
				elif mix_counts.asNumber() == NEG_INF:
					new_media[mol_id] = 0.0 * CONC_UNITS

				else:
					new_counts = base_counts + mix_counts
					new_conc = new_counts / new_volume
					new_media[mol_id] = new_conc

			# if mol_id is not in ingredients, dilute its concentration in new_media
			else:
				base_counts = base_conc * base_media_volume
				new_conc = base_counts / new_volume
				new_media[mol_id] = new_conc

		return new_media

	def make_timeline(self, timeline_str):
		'''
		Make a timeline from a string

		Args:
			timeline_str (str): 'time1 media_id1, time2 media_id2'
		Returns:
			timeline (list[tuple]): a list of tuples with (time (float), media_id (str))

		TODO (Eran) make a parsing expression grammar for this: https://github.com/erikrose/parsimonious
		TODO (Eran) expand capabilities to also pass in ingredients to be added from the prior event
		'''

		timeline = []
		events_str = timeline_str.split(', ')
		for event in events_str:
			time, media = event.split()
			timeline.append((float(time),media))

		return timeline