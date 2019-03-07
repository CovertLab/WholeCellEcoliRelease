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

		self.n_tfs = len(sim_data.process.transcription_regulation.tf_ids)


	# Allocate memory
	def allocate(self):
		super(RnaSynthProb, self).allocate()

		self.rnaSynthProb = np.zeros(len(self.rnaIds), np.float64)
		self.gene_copy_number = np.zeros(len(self.rnaIds), np.int64)

		self.pPromoterBound = np.zeros(self.n_tfs, np.float64)
		self.nPromoterBound = np.zeros(self.n_tfs, np.float64)
		self.nActualBound = np.zeros(self.n_tfs, np.float64)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			rnaIds = list(self.rnaIds),
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
			)
