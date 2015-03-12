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

	def getMass(self, ids):
		assert isinstance(ids, list) or isinstance(ids, np.ndarray)
		idx = [np.where(self._allMass['id'] == re.sub("\[[a-z]\]","", i))[0][0] for i in ids]
		return self._allMass['mass'][idx]

	def _buildAllMasses(self, raw_data, sim_data):
		size = len(raw_data.rnas) + len(raw_data.proteins) + len(raw_data.proteinComplexes) + len(raw_data.metabolites) + len(raw_data.polymerized)
		allMass = np.empty(size,
			dtype = [
					('id',		'a50'),
					('mass',	"f8")
					]
			)

		listMass = []
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.rnas])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.proteins])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.proteinComplexes])
		listMass.extend([(x['id'],np.sum(x['mw7.2'])) for x in raw_data.metabolites])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.polymerized])

		allMass[:] = listMass

		field_units = {
			'id'		:	None,
			'mass'		:	units.g / units.mol,
			}

		self._allMass = UnitStructArray(allMass, field_units)