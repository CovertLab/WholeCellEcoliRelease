"""
SimulationData constants
"""

from __future__ import absolute_import, division, print_function

import scipy.constants

from wholecell.utils import units


class Constants(object):
	""" Constants """

	def __init__(self, raw_data):
		self._build_constants(raw_data)

	def _build_constants(self, raw_data, ):
		self.n_avogadro = scipy.constants.Avogadro / units.mol

		parameters = raw_data.parameters
		self.__dict__.update(parameters)
