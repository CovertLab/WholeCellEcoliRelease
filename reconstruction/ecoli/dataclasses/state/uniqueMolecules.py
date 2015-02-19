"""
SimulationData for unique molecules state

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/19/2015
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class UniqueMolecules(object):
	""" UniqueMolecules """

	def __init__(self, raw_data, sim_data):
		self._buildUniqueMolecules(raw_data, sim_data)


	def _buildUniqueMolecules(self, raw_data, sim_data):
		pas