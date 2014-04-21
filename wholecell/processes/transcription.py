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

	# Constructor
	def __init__(self):
		self.meta = {
		"id": "Transcription",
		"name": "Transcription"
		}
		
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

		self.ntps.requestIs(nPolymerizationReactions // 4)
		self.h2o.requestIs(nPolymerizationReactions)
		self.rnapSubunits.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		rnaPolymerases = (self.rnapSubunits.counts() // [2, 1, 1, 1]).min()

		ntpEstimate = 1.1 * 4 * self.ntps.counts().min()

		enzLimit = np.min([
			ntpEstimate,
			rnaPolymerases * self.elngRate * self.timeStepSec
			])

		newRnas = 0
		ntpsUsed = np.zeros(4)

		nRnasToCreate = int(enzLimit / self.avgRnaLength)

		ntpsShape = self.ntps.counts().shape

		rnasCreated = np.zeros_like(self.rnas.counts())

		while enzLimit > 0 and nRnasToCreate > 0:
			if not np.any(
					np.all(
						self.ntps.counts() > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not np.any(enzLimit > np.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if np.sum(self.rnaSynthProb[np.all(self.ntps.counts() > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if np.sum(self.rnaSynthProb[enzLimit > np.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdxs = np.where(self.randStream.mnrnd(nRnasToCreate, self.rnaSynthProb))[0]

			if np.any(self.ntps.counts() < np.sum(self.rnaNtCounts[newIdxs, :], axis = 0)):
				if nRnasToCreate > 0:
					nRnasToCreate //= 2
					continue
				else:
					break

			if enzLimit < np.sum(np.sum(self.rnaNtCounts[newIdxs, :], axis = 0)):
				if nRnasToCreate > 0:
					nRnasToCreate //= 2
					continue
				else:
					break

			enzLimit -= np.sum(self.rnaNtCounts[newIdxs, :])

			self.ntps.countsDec(
				np.sum(self.rnaNtCounts[newIdxs, :], axis = 0).reshape(ntpsShape)
				)

			self.h2o.countDec(1) # TODO: verify this
			self.ppi.countInc(np.sum(self.rnaLens[newIdxs]))
			self.proton.countInc(1) # TODO: verify this

			rnasCreated[newIdxs] += 1

		self.rnas.countsInc(rnasCreated)