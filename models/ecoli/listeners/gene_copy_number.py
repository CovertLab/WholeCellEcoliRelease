#!/usr/bin/env python

"""
RrnCopyNumber

RrnCopyNumber listener. Tracks gene copy number changing due to replication.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/4/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class RrnCopyNumber(wholecell.listeners.listener.Listener):
	""" RrnCopyNumber """

	_name = 'RrnCopyNumber'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RrnCopyNumber, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RrnCopyNumber, self).initialize(sim, sim_data)

		self.bulkMolecules = sim.states['BulkMolecules']

		self.geneIds = sim_data.process.replication.geneData['name']

		self.geneView = self.bulkMolecules.container.countsView(self.geneIds)

	# Allocate memory
	def allocate(self):
		super(RrnCopyNumber, self).allocate()

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
			simulationStep = self.simulationStep(),
			gene_copy_number = self.gene_copy_number,
			total_copy_number = self.total_copy_number,
			)
