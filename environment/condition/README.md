This directory contains files with data on environmental conditions. There is a corresponding directory ```reconstruction/ecoli/flat/condition/``` that includes files with data on how e. coli responds to some of these conditions.
 
In particular, the media ```minimal```, ```minimal_plus_amino_acids```, ```with_aa```, correspond to data in ```reconstruction/ecoli/flat/condition/```. That corresponding naming convention needs to be maintained for proper model functioning.

# Example use of make_media

### Create media object
> media_obj = Media()

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