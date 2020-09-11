#!/usr/bin/env python

"""
mRNACounts Listener
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import wholecell.listeners.listener


class mRNACounts(wholecell.listeners.listener.Listener):
	"""
	Listener for the counts of each mRNA species. Includes the counts of both
	partial and full transcripts.
	"""
	_name = 'mRNACounts'

	def __init__(self, *args, **kwargs):
		super(mRNACounts, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(mRNACounts, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

		# Get IDs and indexes of all mRNAs
		self.all_RNA_ids = sim_data.process.transcription.rna_data['id']
		self.mRNA_indexes = np.where(
			sim_data.process.transcription.rna_data['is_mRNA'])[0]
		self.mRNA_ids = self.all_RNA_ids[self.mRNA_indexes]

	def allocate(self):
		super(mRNACounts, self).allocate()

		self.mRNA_counts = np.zeros(len(self.mRNA_ids), dtype=np.int64)
		self.full_mRNA_counts = np.zeros(len(self.mRNA_ids), dtype=np.int64)
		self.partial_mRNA_counts = np.zeros(len(self.mRNA_ids), dtype=np.int64)

	def update(self):
		# Get attributes of mRNAs
		RNAs = self.uniqueMolecules.container.objectsInCollection('RNA')
		TU_indexes, can_translate, is_full_transcript = RNAs.attrs(
			'TU_index', 'can_translate', 'is_full_transcript')

		# Get counts of all mRNAs, full mRNAs, and partial mRNAs
		self.mRNA_counts = np.bincount(
			TU_indexes[can_translate],
			minlength=len(self.all_RNA_ids))[self.mRNA_indexes]
		self.full_mRNA_counts = np.bincount(
			TU_indexes[np.logical_and(can_translate, is_full_transcript)],
			minlength=len(self.all_RNA_ids))[self.mRNA_indexes]
		self.partial_mRNA_counts = self.mRNA_counts - self.full_mRNA_counts

	def tableCreate(self, tableWriter):
		subcolumns = {
			'mRNA_counts': 'mRNA_ids',
			'full_mRNA_counts': 'mRNA_ids',
			'partial_mRNA_counts': 'mRNA_ids',
			}

		tableWriter.writeAttributes(
			mRNA_ids = self.mRNA_ids.tolist(),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			mRNA_counts = self.mRNA_counts,
			full_mRNA_counts = self.full_mRNA_counts,
			partial_mRNA_counts = self.partial_mRNA_counts,
			)
