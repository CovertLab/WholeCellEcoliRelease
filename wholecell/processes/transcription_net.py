#!/usr/bin/env python

"""
TranscriptionNet

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/23/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class TranscriptionNet(wholecell.processes.process.Process):
	""" TranscriptionNet """

	_name = "TranscriptionNet"


	# Construct object graph
	def initialize(self, sim, kb):
		super(TranscriptionNet, self).initialize(sim, kb)

		# Views
		self.ntps = self.bulkMoleculesView(["ATP[c]", "UTP[c]", "CTP[c]", "GTP[c]"])
		self.nmps = self.bulkMoleculesView(["AMP[c]", "UMP[c]", "CMP[c]", "GMP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')


	def calculateRequest(self):
		self.ntps.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		counts = self.ntps.counts()

		self.ppi.countInc(counts.sum())

		self.nmps.countsInc(counts)

		self.dntps.countsDec(counts)
