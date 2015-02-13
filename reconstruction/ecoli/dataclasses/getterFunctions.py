"""
SimulationData getter functions

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses.dataclass

import re
import numpy as np

class getterFunctions(reconstruction.ecoli.dataclasses.dataclass.DataClass):
	""" getterFunctions """

	def __init__(self, simData):
		super(getterFunctions, self).__init__(simData)

	def getMass(self, ids):
		assert isinstance(ids, list) or isinstance(ids, np.ndarray)
		idx = [np.where(self._simData._allMass['id'] == re.sub("\[[a-z]\]","", i))[0][0] for i in ids]
		return self._simData._allMass['mass'][idx]