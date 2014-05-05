#!/usr/bin/env python

"""
Replication

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/18/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None

		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		# Load parameters
		dNtpIds = ['DATP[c]', 'DTTP[c]', 'DCTP[c]', 'DGTP[c]']
		dNmpIds = ['DAMP[n]', 'DTMP[n]', 'DCMP[n]', 'DGMP[n]']

		self.sequence = kb.genomeSeq

		# Views
		self.dntps = self.bulkMoleculesView(dNtpIds)
		self.dnmps = self.bulkMoleculesView(dNmpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.dnaPolymerase = self.uniqueMoleculesView('activeDnaPolymerase')

	def calculateRequest(self):
		self.dntps.requestAll()
		self.h2o.requestIs(self.dntps.total().sum())


	# Calculate temporal evolution
	def evolveState(self):
		counts = self.dntps.counts()

		self.ppi.countInc(counts.sum())

		self.dnmps.countsInc(counts)

		self.dntps.countsDec(counts)

		self.h2o.countDec(counts.sum())