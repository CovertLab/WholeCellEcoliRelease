"""
SimulationData for replication process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses.dataclass

class Replication(reconstruction.ecoli.dataclasses.dataclass.DataClass):
	""" Replication """

	def __init__(self, simData):
		super(Replication, self).__init__(simData)
		self._buildSequence()

	def _buildSequence(self):
		self.genome_sequence = self._simData.raw_data.genome_sequence
		self.genome_length = len(self.genome_sequence)
		self.genome_A_count = self.genome_sequence.count("A")
		self.genome_T_count = self.genome_sequence.count("T")
		self.genome_G_count = self.genome_sequence.count("G")
		self.genome_C_count = self.genome_sequence.count("C")