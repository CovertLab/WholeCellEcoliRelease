"""
SimulationData relation functions
"""

from __future__ import absolute_import, division, print_function

import numpy as np


class Relation(object):
	""" Relation """

	def __init__(self, raw_data, sim_data):
		self._build_RNA_to_monomer_mapping(raw_data, sim_data)
		self._build_monomer_to_mRNA_mapping(raw_data, sim_data)

	def _build_RNA_to_monomer_mapping(self, raw_data, sim_data):
		"""
		Build a vector that can map vectors that describe a property for RNAs
		into a vector that describes the same property for the corresponding
		monomers if used as an index array. Assumes that each monomer maps to a
		single RNA (A single RNA can map to multiple monomers).

		e.g.
		monomer_property = RNA_property[
			sim_data.relation.RNA_to_monomer_mapping]
		"""
		# Map RNA IDs to indexes given in rna_data (rnas.tsv)
		rna_id_to_index = {
			rna_id: i for i, rna_id
			in enumerate(sim_data.process.transcription.rna_data['id'])
			}

		# List the rna_data indexes of RNA IDs in the order of corresponding
		# RNAs given in monomer_data (proteins.tsv)
		self.RNA_to_monomer_mapping = np.array([
			rna_id_to_index[rna_id] for rna_id
			in sim_data.process.translation.monomer_data['rna_id']
			])

	def _build_monomer_to_mRNA_mapping(self, raw_data, sim_data):
		"""
		Builds a sparse matrix that can map vectors that describe a property
		for protein monomers into a vector that describes the same property for
		the corresponding mRNAs if multiplied to the right of the original
		vector. The transformed property must be additive (i.e. if two proteins
		map to the same mRNA, the values given for the two proteins are added
		to yield a value for the mRNA).

		The full matrix can be returned by using the monomer_to_mRNA_mapping()
		method.
		"""
		# Initialize sparse matrix variables
		self._monomer_to_mRNA_mapping_i = []
		self._monomer_to_mRNA_mapping_j = []
		self._monomer_to_mRNA_mapping_v = []
		self._monomer_to_mRNA_mapping_shape = (
			len(sim_data.process.translation.monomer_data),
			sim_data.process.transcription.rna_data['is_mRNA'].sum()
		)

		# Build mapping from mRNA ID to mRNA index
		mRNA_data = sim_data.process.transcription.rna_data[
			sim_data.process.transcription.rna_data['is_mRNA']]
		mRNA_id_to_index = {
			mRNA['id']: j for j, mRNA in enumerate(mRNA_data)
		}

		# Build sparse matrix
		for i, monomer in enumerate(sim_data.process.translation.monomer_data):
			self._monomer_to_mRNA_mapping_i.append(i)
			self._monomer_to_mRNA_mapping_j.append(
				mRNA_id_to_index[monomer['rna_id']])
			self._monomer_to_mRNA_mapping_v.append(1)

	def monomer_to_mRNA_mapping(self):
		"""
		Returns the full version of the sparse matrix built by
		_build_monomer_to_mRNA_mapping().

		e.g.
		mRNA_property = monomer_property.dot(
			sim_data.relation.monomer_to_mRNA_mapping())
		"""
		out = np.zeros(self._monomer_to_mRNA_mapping_shape, dtype=np.float64)
		out[self._monomer_to_mRNA_mapping_i, self._monomer_to_mRNA_mapping_j] = self._monomer_to_mRNA_mapping_v
		return out
