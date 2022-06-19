"""
SimulationData for unique molecules state
"""

from __future__ import absolute_import, division, print_function

import collections

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.stateFunctions import addToStateCommon


class UniqueMolecules(object):
	""" UniqueMolecules """

	def __init__(self, raw_data, sim_data):
		self.unique_molecule_definitions = collections.OrderedDict()

		uniqueMoleculeMasses = np.zeros(0,
				dtype = [
						("id", "U50"),
						("mass", "{}f8".format(len(sim_data.submass_name_to_index))),
						]
			)
		field_units = {
			"id" : None,
			"mass" : units.g / units.mol
			}

		self.unique_molecule_masses = UnitStructArray(uniqueMoleculeMasses, field_units)

	def add_to_unique_state(self, uniqueId, attributeDef, mass):
		self.unique_molecule_definitions.update({uniqueId : attributeDef})
		self.unique_molecule_masses = addToStateCommon(self.unique_molecule_masses, [uniqueId], mass)
