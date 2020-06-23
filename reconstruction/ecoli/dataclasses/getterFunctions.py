"""
SimulationData getter functions

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import absolute_import, division, print_function

import itertools
import re
from typing import List, Union

import numpy as np

# Unit imports
from wholecell.utils import units


class getterFunctions(object):
	""" getterFunctions """

	def __init__(self, raw_data, sim_data):
		_ = sim_data  # unused
		self._buildAllMasses(raw_data)
		self._buildLocations(raw_data)

	def getMass(self, ids):
		assert isinstance(ids, (list, np.ndarray))
		masses = [self._all_mass[self._location_tag.sub('', i)] for i in ids]
		return self._mass_units * np.array(masses)

	def getLocation(self, ids):
		# type: (Union[List[str], np.ndarray]) -> List[str]
		assert isinstance(ids, (list, np.ndarray))
		return [self._locationDict[x] for x in ids]

	def get_location_tag(self, id_):
		# type: (str) -> str
		"""Look up a location id and return a location suffix tag like '[c]'."""
		return '[{}]'.format(self._locationDict[id_][0])

	def check_valid_molecule(self, mol_id):
		return mol_id in self._all_mass and mol_id in self._locationDict

	def _buildAllMasses(self, raw_data):
		all_mass = {}
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.rnas})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.proteins})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.proteinComplexes})
		all_mass.update({x['id']: np.sum(x['mw7.2']) for x in raw_data.metabolites})
		all_mass.update({x['id']: np.sum(x['mw7.2']) for x in raw_data.modifiedForms})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.polymerized})
		all_mass.update({x['id']: np.sum(x['mw7.2']) for x in raw_data.water})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.full_chromosome})

		self._all_mass = all_mass
		self._mass_units = units.g / units.mol
		self._location_tag = re.compile(r'\[[a-z]\]')

	def _buildLocations(self, raw_data):
		locationDict = {
			item["id"]: list(item["location"])
			for item in itertools.chain(
				raw_data.rnas,
				raw_data.proteins,
				raw_data.proteinComplexes,
				raw_data.metabolites,
				raw_data.polymerized,
				raw_data.water,
				raw_data.modifiedForms)}

		self._locationDict = locationDict
