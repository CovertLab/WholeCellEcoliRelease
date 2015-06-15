"""
SimulationData relation functions

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/10/2015
"""

from __future__ import division

import re
import numpy as np

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class Relation(object):
	""" Relation """

	def __init__(self, raw_data, sim_data):
		self._buildRnaIndexToMonomerMapping(raw_data, sim_data)
		self._buildMonomerIndexToRnaMapping(raw_data, sim_data)
		#self._buildRnaIndexToGeneMapping(raw_data, sim_data)

	def _buildRnaIndexToMonomerMapping(self, raw_data, sim_data):
		self.rnaIndexToMonomerMapping = np.array([np.where(x == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.translation.monomerData["rnaId"]])

	def _buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		self.monomerIndexToRnaMapping = np.array([np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] for x in sim_data.process.transcription.rnaData["id"] if len(np.where(x == sim_data.process.translation.monomerData["rnaId"])[0])])

	#def _buildRnaIndexToGeneMapping(self, raw_data, sim_data):
	#	self.rnaIndexToGeneMapping = np.array([np.where(x + "[c]" == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.replication.geneData["rnaId"]])