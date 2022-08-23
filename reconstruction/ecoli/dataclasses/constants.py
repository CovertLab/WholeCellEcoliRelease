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
		parameters = raw_data.parameters
		self.__dict__.update(parameters)

		# Add constants not specified in raw_data.parameters
		self.n_avogadro = scipy.constants.Avogadro / units.mol
		self.c_period = len(raw_data.genome_sequence) * units.nt / self.replisome_elongation_rate / 2
