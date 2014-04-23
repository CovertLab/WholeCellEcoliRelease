#!/usr/bin/env python

"""
Transcription

Transcription sub-model. Encodes molecular simulation of macromolecular bacterial transcription

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class Transcription(wholecell.processes.process.Process):
	""" Transcription """

	_name = "Transcription"

	# Constructor
	def __init__(self):
		# Partitions
		self.metabolitePartition = None
		self.rnaPartition = None
		self.enzymePartition = None

		# Constants
		self.cellCycleLength = None
		self.elngRate = None
		self.rnaLens = None					# RNA lengths
		self.rnaNtCounts = None				# RNA nucleotide counts [nt x RNA] <-- TODO: Check this
		self.rnaSynthProb = None			# Relative RNA synthesis rates

		super(Transcription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Transcription, self).initialize(sim, kb)

		# Load parameters
		self.cellCycleLength = kb.cellCycleLen.to('s').magnitude
		self.elngRate = kb.rnaPolymeraseElongationRate.to('nucleotide / s').magnitude

		rnaIds = kb.rnaData['id']

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		# RNA
		self.rnaNtCounts = kb.rnaData['countsAUCG']
		self.rnaLens = kb.rnaData['length']
		self.avgRnaLength = np.mean(self.rnaLens)
		
		self.rnaSynthProb = kb.rnaData['synthProb']

		# Views
		self.ntps = self.bulkMoleculesView(["ATP[c]", "UTP[c]", "CTP[c]", "GTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')

		self.rnas = self.bulkMoleculesView(rnaIds)

		self.rnapSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		rnaPolymerases = (self.rnapSubunits.total() // [2, 1, 1, 1]).min()

		ntpEstimate = 4 * self.ntps.total().min()

		nPolymerizationReactions = np.min([
			ntpEstimate,
			rnaPolymerases * self.elngRate * self.timeStepSec
			])

		self.ntps.requestAll()
		# self.ntps.requestIs(nPolymerizationReactions // 4)
		self.h2o.requestIs(nPolymerizationReactions)
		self.rnapSubunits.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		rnaPolymerases = (self.rnapSubunits.counts() // [2, 1, 1, 1]).min()

		enzLimit = rnaPolymerases * self.elngRate * self.timeStepSec

		from wholecell.utils.simple_polymerize import simplePolymerize

		ntpCounts, rnasCreated = simplePolymerize(
			self.rnaNtCounts,
			enzLimit,
			self.ntps.counts(),
			self.rnaSynthProb,
			self.randStream
			)

		self.ntps.countsIs(ntpCounts)

		self.rnas.countsInc(rnasCreated)

		nRnasCreated = rnasCreated.sum()
		nNtpsUsed = np.dot(rnasCreated, self.rnaNtCounts).sum()

		# This assumes 5' triphosphate is hydrolyzed and 5' OH is deprotonated
		#                      O
		#                     ||
		# End of chain is (-)O-P-O-CH2-Ribose
		#                      |
		#                      O(-)

		self.h2o.countDec(nRnasCreated)
		self.proton.countInc(nRnasCreated)

		self.ppi.countInc(nNtpsUsed)
