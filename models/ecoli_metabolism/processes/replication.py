#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		# Load constants

		dNtpIDs = kb.dNtpIds

		sequence = kb.genomeSeq

		monomerCounts = np.array([sequence.count(s) for s in ("A", "C", "G", "T")], np.float64)

		self.monomerComposition = monomerCounts / monomerCounts.sum()

		self.maxPolymerizationRate = 2 * kb.dnaPolymeraseElongationRate.asNumber() * self.timeStepSec

		self.maxIncorporated = 2 * len(sequence)

		# Create views on state

		self.dntps = self.bulkMoleculesView(dNtpIDs)

		self.ppi = self.bulkMoleculeView("PPI[c]")


	def calculateRequest(self):
		monomersIncoporated = 0 # TODO: constrain elongation by self.maxIncorporated

		monomersRemaining = self.maxIncorporated - monomersIncoporated

		totalMonomers = min(monomersRemaining, self.maxPolymerizationRate)

		dntpsRequested = self.randomState.multinomial(
			totalMonomers,
			self.monomerComposition
			)

		if (dntpsRequested > self.dntps.total()).any():
			# TODO: flag simulation instead of printing
			print "{} is metabolically limited".format(self.name())

		self.dntps.requestIs(dntpsRequested)


	def evolveState(self):
		self.ppi.countInc(self.dntps.counts().sum())
		self.dntps.countsIs(0)
