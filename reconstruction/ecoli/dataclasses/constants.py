"""
SimulationData constants

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/9/2015
"""

from __future__ import absolute_import, division, print_function

import scipy.constants

from wholecell.utils import units


class Constants(object):
	""" Constants """

	def __init__(self, raw_data, sim_data):
		self._buildConstants(raw_data, sim_data)

	def _buildConstants(self, raw_data, sim_data):
		self.nAvogadro = scipy.constants.Avogadro / units.mol

		parameters = raw_data.parameters
		self.__dict__.update(parameters)
