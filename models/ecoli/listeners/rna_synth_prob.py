#!/usr/bin/env python

"""
RnaSynthProb

Records RNA synthesis probabilities

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2016
"""

from __future__ import division

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
		self.rnaIds = sim_data.process.transcription.rnaData["id"]
		self.n_TU = len(self.rnaIds)

		self.tf_ids = sim_data.process.transcription_regulation.tf_ids
		self.n_TF = len(self.tf_ids)


	# Allocate memory
	def allocate(self):
		super(RnaSynthProb, self).allocate()

		self.rnaSynthProb = np.zeros(self.n_TU, np.float64)
		self.gene_copy_number = np.zeros(self.n_TU, np.int16)

		self.pPromoterBound = np.zeros(self.n_TF, np.float64)
		self.nPromoterBound = np.zeros(self.n_TF, np.float64)
		self.nActualBound = np.zeros(self.n_TF, np.float64)

		# This array gets flattened at tableAppend(). Resulting array should
		# be reshaped before use.
		self.n_bound_TF_per_TU = np.zeros((self.n_TU, self.n_TF), np.int16)

		# Properties of bound TFs
		self.bound_TF_indexes = np.array([], np.int64)
		self.bound_TF_coordinates = np.array([], np.int64)
		self.bound_TF_domains = np.array([], np.int64)


	def update(self):
		promoters = self.uniqueMolecules.container.objectsInCollection('promoter')

		if len(promoters) > 0:
			TU_indexes, all_coordinates, all_domains, bound_TFs = promoters.attrs(
				"TU_index", "coordinates", "domain_index", "bound_TF"
				)

			self.gene_copy_number = np.bincount(TU_indexes, minlength=self.n_TU)

			bound_promoter_indexes, TF_indexes = np.where(bound_TFs)

			self.bound_TF_indexes = TF_indexes
			self.bound_TF_coordinates = all_coordinates[bound_promoter_indexes]
			self.bound_TF_domains = all_domains[bound_promoter_indexes]
		else:
			self.bound_TF_indexes = np.array([])
			self.bound_TF_coordinates = np.array([])
			self.bound_TF_domains = np.array([])


	def tableCreate(self, tableWriter):
		subcolumns = {
			'gene_copy_number': 'rnaIds',
			'rnaSynthProb': 'rnaIds',
			'pPromoterBound': 'tf_ids',
			'nPromoterBound': 'tf_ids',
			'nActualBound': 'tf_ids',
			'n_bound_TF_per_TU': 'rnaIds'}

		tableWriter.writeAttributes(
			rnaIds = list(self.rnaIds),
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
			gene_copy_number = self.gene_copy_number,
			pPromoterBound = self.pPromoterBound,
			nPromoterBound = self.nPromoterBound,
			nActualBound = self.nActualBound,
			n_bound_TF_per_TU = self.n_bound_TF_per_TU,
			bound_TF_indexes = self.bound_TF_indexes,
			bound_TF_coordinates = self.bound_TF_coordinates,
			bound_TF_domains = self.bound_TF_domains,
			)
