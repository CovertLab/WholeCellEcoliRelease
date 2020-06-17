'''
Functions for making media

# Example use of make_media

### Create media object
> media_obj = Media(raw_data)

### Retrieve stock media
> base_media = media_obj.stock_media['M9_GLC']
> base_media2 = media_obj.stock_media['5X_supplement_EZ']

### Define a dict of ingredients
Ingredients is a dict with molecule ids as the keys.
Each ingredient's value is a dict with {'weight': value * (units.g), 'counts': value * (units.mmol), 'volume': value *  (units.L)}.
Only one of 'weights' (in units.g) or 'counts' (in units.mmol) is required; if both are specified, it will use weight.
If weight or counts is Infinity, it sets the final concentration to inf.
If weight or counts is -Infinity, it sets the final concentration to 0.

Example:
> ingredients = {
	'L-ALPHA-ALANINE': {'weight': 1.78 * units.g, 'volume': 0.025 * units.L},
	'ARG': {'weight': 8.44 * units.g, 'volume': 0.1 * units.L},
	'UREA': {'counts': 102.0 * units.mmol, 'volume': 1.0 * units.L},
	'LEU': {'weight': float("inf") * units.g, 'volume': 0 * units.L},
	'OXYGEN-MOLECULE': {'weight': float("-inf") * units.g, 'volume': 0 * units.L},
    }

### Add ingredients directly into an existing media
> new_media1 = media_obj.add_ingredients(base_media, 0.8 * units.L, ingredients)

### Combine two medias
> new_media2 = media_obj.combine_media(base_media, 0.8 * units.L, base_media2, 0.2 * units.L)

'''

from __future__ import absolute_import, division, print_function

from wholecell.utils import units
import six


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

	def __init__(self, raw_data):
		# get dicts from knowledge base
		self.environment_molecules_fw = self._get_environment_molecules_fw(raw_data)
		self.stock_media = self._get_stock_media(raw_data)
		self.recipes = self._get_recipes(raw_data)

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

	def _get_recipes(self, raw_data):
		'''load recipes'''

		recipes = {}
		for row in raw_data.condition.media_recipes:
			new_media_id = row["media id"]

			recipe = {}
			recipe["base media"] = row["base media"]
			recipe["added media"] = row.get("added media", None)
			recipe["ingredients"] = row.get("ingredients", None)

			recipe["base media volume"] = row.get("base media volume", 0 * units.L)
			recipe["added media volume"] = row.get("added media volume", 0 * units.L)
			recipe["ingredients weight"] = row.get("ingredients weight", None)
			recipe["ingredients counts"] = row.get("ingredients counts", None)
			recipe["ingredients volume"] = row.get("ingredients volume", 0 * units.L)

			recipes[new_media_id] = recipe

		return recipes

	def make_saved_media(self):
		'''make all the media recipes in self.recipes'''

		self.saved_media = {}
		for new_media_id in self.recipes:
			new_media = self.make_recipe(new_media_id)
			self.saved_media[new_media_id] = new_media

		return self.saved_media

	def make_recipe(self, media_id):
		'''make a single media recipe from self.recipes'''

		recipe = self.recipes[media_id]
		base_id = recipe["base media"]
		added_media_id = recipe["added media"]
		ingredient_ids = recipe["ingredients"]
		base_media = self.stock_media[base_id]
		base_vol = recipe["base media volume"]

		if added_media_id:
			added_media = self.stock_media[added_media_id]
			added_vol = recipe["added media volume"]
			new_media = self.combine_media(base_media, base_vol, added_media, added_vol)
			base_media = new_media
			base_vol += added_vol

		if ingredient_ids:
			added_weight = recipe.get("ingredients weight", [])
			added_counts = recipe.get("ingredients counts", [])
			added_vol = recipe.get("ingredients volume")  # the row is a list with units.L, even an empty list is read.
			ingredients = {ingred_id: {} for ingred_id in ingredient_ids}
			for index, ingred_id in enumerate(ingredient_ids):
				if len(added_weight):
					ingredients[ingred_id]['weight'] = added_weight[index]
				if len(added_counts):
					ingredients[ingred_id]['counts'] = added_counts[index]
				if len(added_vol):
					ingredients[ingred_id]['volume'] = added_vol[index]
				else:
					ingredients[ingred_id]['volume'] = 0 * units.L
			new_media = self.add_ingredients(base_media, base_vol, ingredients)

		if not added_media_id and not ingredient_ids:
			new_media = base_media

		# remove concentration units, setting at CONC_UNITS
		unitless_new_media = {mol: conc.asNumber(CONC_UNITS) for mol, conc in six.viewitems(new_media)}

		return unitless_new_media

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
		new_media = {mol_id: 0.0 * CONC_UNITS for mol_id, conc in six.viewitems(base_media)}

		# get new_media volume
		new_volume = base_media_volume + mix_media_volume

		for mol_id, base_conc in six.viewitems(base_media):
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
		new_media = {mol_id: 0.0 * CONC_UNITS for mol_id, conc in six.viewitems(base_media)}

		# get new_media volume
		ingredients_volume = 0 * VOLUME_UNITS
		for quantities in six.viewvalues(ingredients):
			ingredients_volume += quantities['volume']
		new_volume = base_media_volume + ingredients_volume

		# get new_media concentrations from mixing ingredients
		for mol_id, base_conc in six.viewitems(base_media):

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
