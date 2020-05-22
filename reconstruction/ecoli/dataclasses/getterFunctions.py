"""
SimulationData getter functions

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import re
import numpy as np

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class getterFunctions(object):
	""" getterFunctions """

	def __init__(self, raw_data, sim_data):
		self._buildAllMasses(raw_data, sim_data)
		self._buildLocations(raw_data, sim_data)

	def getMass(self, ids):
		assert isinstance(ids, (list, np.ndarray))
		try:
			masses = [self._all_mass[self._location_tag.sub('', i)] for i in ids]
		except KeyError:
			raise Exception("Unrecognized id: {}".format(i))

		return self._mass_units * np.array(masses)

	def getLocation(self, ids):
		assert isinstance(ids, list) or isinstance(ids, np.ndarray)
		return [self._locationDict[x] for x in ids]

	def check_valid_molecule(self, mol_id):
		return mol_id in self._all_mass and mol_id in self._locationDict

	def _buildAllMasses(self, raw_data, sim_data):
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
		self._location_tag = re.compile('\[[a-z]\]')

	def _buildLocations(self, raw_data, sim_data):
		locationDict = {}
		for item in raw_data.rnas:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]
		for item in raw_data.proteins:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]
		for item in raw_data.proteinComplexes:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]
		for item in raw_data.metabolites:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]
		for item in raw_data.polymerized:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]
		for item in raw_data.water:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]
		for item in raw_data.modifiedForms:
			locationDict[item["id"]] = [x.encode("utf-8") for x in item["location"]]

		self._locationDict = locationDict
