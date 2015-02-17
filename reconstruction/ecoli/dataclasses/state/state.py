"""
SimulationData state associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses.dataclass
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkMolecules import BulkMolecules
from reconstruction.ecoli.dataclasses.state.bulkChromosome import BulkChromosome


import re
import numpy as np

class State(reconstruction.ecoli.dataclasses.dataclass.DataClass):
	""" State """

	def __init__(self, simData):
		super(State, self).__init__(simData)

		self.bulkMolecules = BulkMolecules(simData)
		self.bulkChromosome = BulkChromosome(simData)
		self._buildCompartments()


	def _buildCompartments(self):
		compartmentData = np.empty(len(self._simData.raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in self._simData.raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in self._simData.raw_data.compartments]
		self.compartments = compartmentData
