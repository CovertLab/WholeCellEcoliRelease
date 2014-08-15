#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

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
		polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedDNT_IDs]

		sequence = kb.genomeSeq

		monomerCounts = np.array([sequence.count(s) for s in ("A", "C", "G", "T")], np.float64)

		self.monomerComposition = monomerCounts / monomerCounts.sum()

		self.maxPolymerizationRate = 2 * kb.dnaPolymeraseElongationRate.magnitude * self.timeStepSec

		self.maxIncorporated = 2 * len(sequence)

		# Create views on state

		self.polymerized = self.bulkMoleculesView(polymerizedIDs)

		self.dntps = self.bulkMoleculesView(dNtpIDs)

		self.ppi = self.bulkMoleculeView("PPI[c]")


	def calculateRequest(self):
		monomersIncoporated = self.polymerized.total().sum()

		monomersRemaining = max(self.maxIncorporated - monomersIncoporated, 0)

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
		dntpCounts = self.dntps.counts()

		self.polymerized.countsInc(dntpCounts)

		self.ppi.countInc(dntpCounts.sum())
		self.dntps.countsIs(0)

		print dntpCounts.sum()
