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
import tables

import wholecell.listeners.listener
from reconstruction.ecoli.fitter import normalize

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

		self.geneIds = kb.geneData['name']

		self.geneView = self.bulkChromosome.container.countsView(self.geneIds)

	# Allocate memory
	def allocate(self):
		super(GeneCopyNumber, self).allocate()

		self.gene_copy_number = np.zeros(len(self.geneIds))
		self.total_copy_number = 0

	def update(self):
		self.gene_copy_number = self.geneView.counts()
		self.total_copy_number = np.sum(self.gene_copy_number)

	def pytablesCreate(self, h5file, expectedRows):

		# Columns
		d = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"gene_copy_number": tables.UInt64Col(self.gene_copy_number.shape),
			"total_copy_number": tables.UInt64Col(),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self._name,
			d,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		# Store units as metadata
		t.attrs.gene_copy_number = self.countUnits
		t.attrs.total_copy_number = self.countUnits


	def pytablesAppend(self, h5file):

		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["gene_copy_number"] = self.gene_copy_number
		entry["total_copy_number"] = self.total_copy_number

		entry.append()

		t.flush()