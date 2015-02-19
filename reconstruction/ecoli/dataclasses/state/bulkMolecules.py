"""
SimulationData for bulk molecules state

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class BulkMolecules(object):
	""" BulkMolecules """

	def __init__(self, raw_data, sim_data):
		self.bulkMoleculesData = np.zeros(
			0,
			dtype = [
				("id", "a50"),
				("mass", "{}f8".format(len(sim_data.molecular_weight_order))),
				]
			)