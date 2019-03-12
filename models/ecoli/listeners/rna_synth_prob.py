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
		self.n_trs_units = len(self.rnaIds)

		self.tf_ids = sim_data.process.transcription_regulation.tf_ids
		self.n_tfs = len(self.tf_ids)


	# Allocate memory
	def allocate(self):
		super(RnaSynthProb, self).allocate()

		self.rnaSynthProb = np.zeros(self.n_trs_units, np.float64)
		self.gene_copy_number = np.zeros(self.n_trs_units, np.int16)

		self.pPromoterBound = np.zeros(self.n_tfs, np.float64)
		self.nPromoterBound = np.zeros(self.n_tfs, np.float64)
		self.nActualBound = np.zeros(self.n_tfs, np.float64)

		# This array gets flattened at tableAppend(). Resulting array should
		# be reshaped before use.
		self.n_bound_tfs_per_trs_unit = np.zeros(
			(self.n_trs_units, self.n_tfs), np.int16
			)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			rnaIds = list(self.rnaIds),
			tf_ids = list(self.tf_ids),
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
			n_bound_tfs_per_trs_unit = self.n_bound_tfs_per_trs_unit,
			)
