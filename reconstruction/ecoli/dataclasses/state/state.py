"""
SimulationData state associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkMolecules import BulkMolecules
from reconstruction.ecoli.dataclasses.state.bulkChromosome import BulkChromosome

import re
import numpy as np

class State(object):
	""" State """

	def __init__(self, raw_data, sim_data):

		self.bulkMolecules = BulkMolecules(raw_data, sim_data)
		self.bulkChromosome = BulkChromosome(raw_data, sim_data)
		self._buildCompartments(raw_data, sim_data)


	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData
