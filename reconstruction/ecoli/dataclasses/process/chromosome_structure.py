"""
SimulationData for the ChromosomeStructure process

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/10/2020
"""

from __future__ import division, print_function, absolute_import


class ChromosomeStructure(object):
	"""
	SimulationData for the ChromosomeStructure process
	"""

	def __init__(self, raw_data, sim_data):
		self._build_supercoiling_parameters(raw_data, sim_data)


	def _build_supercoiling_parameters(self, raw_data, sim_data):
		"""
		Load parameters used for DNA supercoiling from raw_data.
		"""
		for parameter in raw_data.dna_supercoiling:
			self.__setattr__(parameter['name'], parameter['value'])

		# Unique index used for dummy molecule at terC
		self.terC_dummy_molecule_index = -1
