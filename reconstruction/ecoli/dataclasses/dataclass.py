"""
SimulationData dataclass base class

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

class DataClass(object):
	""" DataClass """

	def __init__(self, simData):
		self._simData = simData