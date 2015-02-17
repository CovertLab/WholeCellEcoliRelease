"""
SimulationData getter functions

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import re
import numpy as np

class getterFunctions(object):
	""" getterFunctions """

	def __init__(self, raw_data, sim_data):
		pass

	def getMass(self, ids):
		assert isinstance(ids, list) or isinstance(ids, np.ndarray)
		idx = [np.where(self._simData._allMass['id'] == re.sub("\[[a-z]\]","", i))[0][0] for i in ids]
		return self._simData._allMass['mass'][idx]