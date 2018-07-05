"""
Simulation data for Environment state

The environmental state object

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.stateFunctions import addToStateCommon

class Environment(object):
	""" Environment """

	def __init__(self, raw_data, sim_data):
		environmentData = np.zeros(
			0,
			dtype = [
				("id", "a50"),
				]
			)

		# Add units to values
		field_units = {
			"id"		:	None,
			}

		self.environmentData = UnitStructArray(environmentData, field_units)
