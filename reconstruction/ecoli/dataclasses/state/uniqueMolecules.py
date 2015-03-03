"""
SimulationData for unique molecules state

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/19/2015
"""

from __future__ import division
import collections

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkStateFunctions import addToStateCommon


class UniqueMolecules(object):
	""" UniqueMolecules """

	def __init__(self, raw_data, sim_data):
		self.uniqueMoleculeDefinitions = collections.OrderedDict()

		uniqueMoleculeMasses = np.zeros(0,
				dtype = [
						("id", "a50"),
						("mass", "{}f8".format(len(sim_data.molecular_weight_order))),
						]
			)
		field_units = {
			"id" : None,
			"mass" : units.g / units.mol
			}

		self.uniqueMoleculeMasses = UnitStructArray(uniqueMoleculeMasses, field_units)

	def addToUniqueState(self, uniqueId, attributeDef, mass):
		self.uniqueMoleculeDefinitions.update({uniqueId : attributeDef})

		self.uniqueMoleculeMasses = addToStateCommon(self.uniqueMoleculeMasses, [uniqueId], mass)