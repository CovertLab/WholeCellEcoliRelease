"""
RnaSynthProb

Records RNA synthesis probabilities
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener


class RnaSynthProb(wholecell.listeners.listener.Listener):
	""" RnaSynthProb """

	_name = "RnaSynthProb"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnaSynthProb, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaSynthProb, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

		self.transcriptInitiation = sim.processes["TranscriptInitiation"]
		self.rna_ids = sim_data.process.transcription.rna_data["id"]
		self.n_TU = len(self.rna_ids)
		self.cistron_ids = sim_data.process.transcription.cistron_data["id"]
		self.n_cistron = len(self.cistron_ids)
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

		self.tf_ids = sim_data.process.transcription_regulation.tf_ids
		self.n_TF = len(self.tf_ids)


	# Allocate memory
	def allocate(self):
		super(RnaSynthProb, self).allocate()

		self.rnaSynthProb = np.zeros(self.n_TU, np.float64)
		self.gene_copy_number = np.zeros(self.n_TU, np.int16)
		self.rna_synth_prob_per_cistron = np.zeros(self.n_cistron, np.float64)

		self.pPromoterBound = np.zeros(self.n_TF, np.float64)
		self.nPromoterBound = np.zeros(self.n_TF, np.float64)
		self.nActualBound = np.zeros(self.n_TF, np.float64)
		self.n_available_promoters = np.zeros(self.n_TF, np.float64)

		# This array gets flattened at tableAppend(). Resulting array should
		# be reshaped before use.
		self.n_bound_TF_per_TU = np.zeros((self.n_TU, self.n_TF), np.int16)

		# Properties of bound TFs
		self.bound_TF_indexes = np.array([], np.int64)
		self.bound_TF_coordinates = np.array([], np.int64)
		self.bound_TF_domains = np.array([], np.int64)


	def update(self):
		promoters = self.uniqueMolecules.container.objectsInCollection('promoter')
		TU_indexes, all_coordinates, all_domains, bound_TFs = promoters.attrs(
			"TU_index", "coordinates", "domain_index", "bound_TF"
			)

		self.gene_copy_number = np.bincount(TU_indexes, minlength=self.n_TU)

		bound_promoter_indexes, TF_indexes = np.where(bound_TFs)

		self.bound_TF_indexes = TF_indexes
		self.bound_TF_coordinates = all_coordinates[bound_promoter_indexes]
		self.bound_TF_domains = all_domains[bound_promoter_indexes]

		self.rna_synth_prob_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.rnaSynthProb)


	def tableCreate(self, tableWriter):
		subcolumns = {
			'gene_copy_number': 'rnaIds',
			'rnaSynthProb': 'rnaIds',
			'rna_synth_prob_per_cistron': 'cistron_ids',
			'pPromoterBound': 'tf_ids',
			'nPromoterBound': 'tf_ids',
			'nActualBound': 'tf_ids',
			'n_available_promoters': 'tf_ids',
			'n_bound_TF_per_TU': 'rnaIds'}

		tableWriter.writeAttributes(
			rnaIds = list(self.rna_ids),
			cistron_ids = list(self.cistron_ids),
			tf_ids = list(self.tf_ids),
			subcolumns = subcolumns)

		tableWriter.set_variable_length_columns(
			'bound_TF_indexes',
			'bound_TF_coordinates',
			'bound_TF_domains',
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			rnaSynthProb = self.rnaSynthProb,
			rna_synth_prob_per_cistron = self.rna_synth_prob_per_cistron,
			gene_copy_number = self.gene_copy_number,
			pPromoterBound = self.pPromoterBound,
			nPromoterBound = self.nPromoterBound,
			nActualBound = self.nActualBound,
			n_available_promoters = self.n_available_promoters,
			n_bound_TF_per_TU = self.n_bound_TF_per_TU,
			bound_TF_indexes = self.bound_TF_indexes,
			bound_TF_coordinates = self.bound_TF_coordinates,
			bound_TF_domains = self.bound_TF_domains,
			)
