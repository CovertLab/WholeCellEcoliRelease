#!/usr/bin/env python

"""
RnaSynthProb

Records RNA synthesis probabilities

@author: Derek Macklin
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
			)
