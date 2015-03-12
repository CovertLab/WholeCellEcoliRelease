"""
SimulationData for replication process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import numpy as np

class Replication(object):
	""" Replication """

	def __init__(self, raw_data, sim_data):
		self._buildSequence(raw_data, sim_data)
		self._buildGeneData(raw_data, sim_data)

	def _buildSequence(self, raw_data, sim_data):
		self.genome_sequence = raw_data.genome_sequence
		self.genome_length = len(self.genome_sequence)
		self.genome_A_count = self.genome_sequence.count("A")
		self.genome_T_count = self.genome_sequence.count("T")
		self.genome_G_count = self.genome_sequence.count("G")
		self.genome_C_count = self.genome_sequence.count("C")

	def _buildGeneData(self, raw_data, sim_data):
		genomeLength = len(raw_data.genome_sequence)
		self.geneData = np.zeros(len(raw_data.genes),
			dtype = [('name'				,	'a50'),
					('rnaId'                ,   'a50'),
					('endCoordinate'		,	'int64')])

		self.geneData['name'] = [x['id'] for x in raw_data.genes]
		self.geneData['rnaId'] = [x['rnaId'] for x in raw_data.genes]
		self.geneData['endCoordinate'] = [(x['coordinate'] + x['length']) % genomeLength if x['direction'] == '+' else (x['coordinate'] - x['length']) % genomeLength for x in raw_data.genes]
