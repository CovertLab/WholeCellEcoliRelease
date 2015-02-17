"""
SimulationData process associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.process.replication import Replication

import re
import numpy as np

class Process(object):
	""" Process """

	def __init__(self, raw_data, sim_data):

		self.replication = Replication(raw_data, sim_data)
