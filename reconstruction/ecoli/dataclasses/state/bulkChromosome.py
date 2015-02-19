"""
SimulationData for bulk chromosome state

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
from reconstruction.ecoli.dataclasses.state.bulkStateFunctions import addToBulkState

class BulkChromosome(object):
	""" BulkChromosome """

	def __init__(self, raw_data, sim_data):
		self.addToBulkState = addToBulkState
		self._buildBulkChromosome(raw_data, sim_data)

	def _buildBulkChromosome(self, raw_data, sim_data):
		bulkChromosome = np.zeros(0,
			dtype = [("id", 			"a50"),
					("mass", "{}f8".format(len(sim_data.molecular_weight_order)))
					]
					)

		# Set genes
		geneIds = [x['id'] for x in raw_data.genes]
		geneMasses = np.zeros((len(geneIds), len(sim_data.molecular_weight_order)), np.float64)

		bulkChromosome = self.addToBulkState(bulkChromosome, geneIds, geneMasses)

		# Add units to values
		field_units = {
			"id"			:	None,
			"mass"					:	units.g / units.mol,
			}
		self.bulkChromosome = UnitStructArray(bulkChromosome, field_units)