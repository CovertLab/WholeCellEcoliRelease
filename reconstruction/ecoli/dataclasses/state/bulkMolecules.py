"""
SimulationData for bulk molecules state
"""

from __future__ import absolute_import, division, print_function

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.stateFunctions import addToStateCommon

class BulkMolecules(object):
	""" BulkMolecules """

	def __init__(self, raw_data, sim_data):
		bulkData = np.zeros(
			0,
			dtype = [
				("id", "U50"),
				("mass", "{}f8".format(len(sim_data.submass_name_to_index))),
				]
			)

		# Add units to values
		field_units = {
			"id"		:	None,
			"mass"				:	units.g / units.mol,
			}

		self.bulk_data = UnitStructArray(bulkData, field_units)

	def add_to_bulk_state(self, ids, masses):
		self.bulk_data = addToStateCommon(self.bulk_data, ids, masses)