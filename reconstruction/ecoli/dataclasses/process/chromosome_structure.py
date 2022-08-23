"""
SimulationData for the ChromosomeStructure process
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
		for name, value in raw_data.dna_supercoiling.items():
			self.__setattr__(name, value)

		# Unique index used for dummy molecule at terC
		self.terC_dummy_molecule_index = -1
