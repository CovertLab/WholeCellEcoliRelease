#!/usr/bin/env python

"""
mRNACounts Listener
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import wholecell.listeners.listener


class mRNACounts(wholecell.listeners.listener.Listener):
	"""
	Listener for the counts of each mRNA transcription units and cistrons.
	Includes the counts of both partial and full transcripts.
	"""
	_name = 'mRNACounts'

	def __init__(self, *args, **kwargs):
		super(mRNACounts, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(mRNACounts, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

		# Get IDs and indexes of all mRNA transcription units
		self.all_TU_ids = sim_data.process.transcription.rna_data['id']
		self.mRNA_indexes = np.where(
			sim_data.process.transcription.rna_data['is_mRNA'])[0]
		self.mRNA_TU_ids = self.all_TU_ids[self.mRNA_indexes]

		# Get IDs and indexes of all mRNA cistrons
		self.all_cistron_ids = sim_data.process.transcription.cistron_data['id']
		self.cistron_is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		self.mRNA_cistron_ids = self.all_cistron_ids[self.cistron_is_mRNA]

		# Get mapping matrix between TUs and cistrons
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

	def allocate(self):
		super(mRNACounts, self).allocate()

		self.mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.full_mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.partial_mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.mRNA_cistron_counts = np.zeros(len(self.mRNA_cistron_ids), dtype=np.int64)

	def update(self):
		# Get attributes of mRNAs
		RNAs = self.uniqueMolecules.container.objectsInCollection('RNA')
		TU_indexes, can_translate, is_full_transcript = RNAs.attrs(
			'TU_index', 'can_translate', 'is_full_transcript')

		# Get counts of mRNA transcription units
		all_TU_counts = np.bincount(
			TU_indexes[can_translate], minlength=len(self.all_TU_ids))
		self.mRNA_counts = all_TU_counts[self.mRNA_indexes]
		full_TU_counts = np.bincount(
			TU_indexes[np.logical_and(can_translate, is_full_transcript)],
			minlength=len(self.all_TU_ids))
		self.full_mRNA_counts = full_TU_counts[self.mRNA_indexes]
		partial_TU_counts = all_TU_counts - full_TU_counts
		self.partial_mRNA_counts = self.mRNA_counts - self.full_mRNA_counts

		# Calculate counts of mRNA cistrons from transcription unit counts
		# TODO (ggsun): Partial mRNA cistron counts should take into account
		# 	the lengths of each mRNA transcript.
		self.mRNA_cistron_counts = self.cistron_tu_mapping_matrix.dot(
			all_TU_counts)[self.cistron_is_mRNA]
		self.full_mRNA_cistron_counts = self.cistron_tu_mapping_matrix.dot(
			full_TU_counts)[self.cistron_is_mRNA]
		self.partial_mRNA_cistron_counts = self.cistron_tu_mapping_matrix.dot(
			partial_TU_counts)[self.cistron_is_mRNA]

	def tableCreate(self, tableWriter):
		subcolumns = {
			'mRNA_counts': 'mRNA_ids',
			'full_mRNA_counts': 'mRNA_ids',
			'partial_mRNA_counts': 'mRNA_ids',
			}

		tableWriter.writeAttributes(
			mRNA_ids = self.mRNA_TU_ids.tolist(),
			mRNA_cistron_ids = self.mRNA_cistron_ids.tolist(),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			mRNA_counts = self.mRNA_counts,
			full_mRNA_counts = self.full_mRNA_counts,
			partial_mRNA_counts = self.partial_mRNA_counts,
			mRNA_cistron_counts = self.mRNA_cistron_counts,
			full_mRNA_cistron_counts = self.full_mRNA_cistron_counts,
			partial_mRNA_cistron_counts = self.partial_mRNA_cistron_counts,
			)
