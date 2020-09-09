"""
SimulationData relation functions
"""

from __future__ import absolute_import, division, print_function

import numpy as np


class Relation(object):
	""" Relation """

	def __init__(self, raw_data, sim_data):
		self._build_rna_index_to_monomer_mapping(raw_data, sim_data)
		self._build_monomer_index_to_rna_mapping(raw_data, sim_data)
		#self._buildRnaIndexToGeneMapping(raw_data, sim_data)

	def _build_rna_index_to_monomer_mapping(self, raw_data, sim_data):
		self.rna_index_to_monomer_mapping = np.array([np.where(x == sim_data.process.transcription.rna_data["id"])[0][0] for x in sim_data.process.translation.monomer_data['rna_id']])

	def _build_monomer_index_to_rna_mapping(self, raw_data, sim_data):
		self.monomer_index_to_rna_mapping = np.array([np.where(x == sim_data.process.translation.monomer_data['rna_id'])[0][0] for x in sim_data.process.transcription.rna_data["id"] if len(np.where(x == sim_data.process.translation.monomer_data['rna_id'])[0])])
