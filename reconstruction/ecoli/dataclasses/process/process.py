"""
SimulationData process associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses.dataclass
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.process.replication import Replication

import re
import numpy as np

class Process(reconstruction.ecoli.dataclasses.dataclass.DataClass):
	""" Process """

	def __init__(self, simData):
		super(Process, self).__init__(simData)

		self.replication = Replication(self._simData)
