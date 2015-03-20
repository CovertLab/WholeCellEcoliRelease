#!/usr/bin/env python

"""
GeneCopyNumber

GeneCopyNumber listener. Tracks gene copy number changing due to replication.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/4/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class GeneCopyNumber(wholecell.listeners.listener.Listener):
	""" GeneCopyNumber """

	_name = 'GeneCopyNumber'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(GeneCopyNumber, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, kb):
		super(GeneCopyNumber, self).initialize(sim, kb)

		self.bulkChromosome = sim.states['BulkChromosome']

		self.geneIds = kb.process.replication.geneData['name']

		self.geneView = self.bulkChromosome.container.countsView(self.geneIds)

	# Allocate memory
	def allocate(self):
		super(GeneCopyNumber, self).allocate()

		self.gene_copy_number = np.zeros(len(self.geneIds))
		self.total_copy_number = 0

	def update(self):
		self.gene_copy_number = self.geneView.counts()
		self.total_copy_number = np.sum(self.gene_copy_number)

	def tableCreate(self, tableWriter):
		# Store units as metadata
		tableWriter.writeAttributes( # TODO: reconsider attribute names
			gene_copy_number = self.countUnits,
			total_copy_number = self.countUnits
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			gene_copy_number = self.gene_copy_number,
			total_copy_number = self.total_copy_number,
			)
