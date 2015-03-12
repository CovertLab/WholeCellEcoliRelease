"""
SimulationData constants

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/9/2015
"""

from __future__ import division


class Constants(object):
	""" Constants """

	def __init__(self, raw_data, sim_data):
		self._buildConstants(raw_data, sim_data)

	def _buildConstants(self, raw_data, sim_data):
		constants = raw_data.constants
		self.__dict__.update(constants)