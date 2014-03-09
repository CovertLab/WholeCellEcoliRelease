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
		self.cellCycleLength = 1 * 3600.	# s # TOKB
		self.elngRate = 50					# nt/s # TOKB
		self.rnaLens = None					# RNA lengths
		self.rnaNtCounts = None				# RNA nucleotide counts [nt x RNA] <-- TODO: Check this
		self.rnaSynthProb = None			# Relative RNA synthesis rates

		super(Transcription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Transcription, self).initialize(sim, kb)

		mc = sim.states["BulkMolecules"]

		rnaIds = [x["id"] + "[c]" for x in kb.rnas]

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.bulkMoleculesPartition.initialize(_metIds + rnaIds + enzIds)

		# Metabolites
		self.bulkMoleculesPartition.metabolites = self.bulkMoleculesPartition.countsView(_metIds)

		self.ntpView = mc.countsView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.bulkMoleculesPartition.ntps = self.bulkMoleculesPartition.countsView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.bulkMoleculesPartition.ppiMol = self.bulkMoleculesPartition.countView('PPI[c]')
		self.bulkMoleculesPartition.h2oMol = self.bulkMoleculesPartition.countView('H2O[c]')
		self.bulkMoleculesPartition.hMol = self.bulkMoleculesPartition.countView('H[c]')

		# RNA
		self.bulkMoleculesPartition.rnas = self.bulkMoleculesPartition.countsView(rnaIds)
		self.rnaNtCounts = np.array([x["ntCount"] for x in kb.rnas])
		self.rnaLens = np.sum(self.rnaNtCounts, axis = 1)
		
		halflives = np.array([x["halfLife"] for x in kb.rnas])
		self.rnaSynthProb = mc._rnaExp * (np.log(2) / self.cellCycleLength + 1 / halflives)
		self.rnaSynthProb /= np.sum(self.rnaSynthProb)

		# Enzymes
		# self.enzyme = sim.states["BulkMolecules"].addPartition(self, ["RNAP70-CPLX[c]"], self.calcReqEnzyme)
		# self.enzyme.idx["rnaPol"] = self.enzyme.getIndex(["RNAP70-CPLX[c]"])[0]
		self.bulkMoleculesPartition.enzymes = self.bulkMoleculesPartition.countsView(enzIds)
		self.bulkMoleculesPartition.rpoAMol = self.bulkMoleculesPartition.countView('EG10893-MONOMER[c]')
		self.bulkMoleculesPartition.rpoBMol = self.bulkMoleculesPartition.countView('RPOB-MONOMER[c]')
		self.bulkMoleculesPartition.rpoCMol = self.bulkMoleculesPartition.countView('RPOC-MONOMER[c]')
		self.bulkMoleculesPartition.rpoDMol = self.bulkMoleculesPartition.countView('RPOD-MONOMER[c]')

		self.rpoAMol = mc.countView('EG10893-MONOMER[c]')
		self.rpoBMol = mc.countView('RPOB-MONOMER[c]')
		self.rpoCMol = mc.countView('RPOC-MONOMER[c]')
		self.rpoDMol = mc.countView('RPOD-MONOMER[c]')

		# Views
		# self.metabolites = self.bulkMoleculesView(_metIds)
		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.rnapSubunits = self.bulkMoleculesView(enzIds)


	def requestBulkMolecules(self):
		self.bulkMoleculesPartition.ntps.countsIs(
			np.min([
				calcRnaps(
					self.rpoAMol.count(), self.rpoBMol.count(),
					self.rpoCMol.count(), self.rpoDMol.count()
					) * self.elngRate * self.timeStepSec,
				4 * np.min(self.ntpView.counts())
				])/4
			)

		self.bulkMoleculesPartition.h2oMol.countIs(1)

		self.bulkMoleculesPartition.rnas.countsIs(0)

		self.bulkMoleculesPartition.enzymes.countsIs(1)


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
		enzLimit = np.min([
			calcRnaps(
				self.bulkMoleculesPartition.rpoAMol.count(),
				self.bulkMoleculesPartition.rpoBMol.count(),
				self.bulkMoleculesPartition.rpoCMol.count(),
				self.bulkMoleculesPartition.rpoDMol.count()
				) * self.elngRate * self.timeStepSec,
			1.1 * 4 * np.min(self.bulkMoleculesPartition.ntps.counts())
			])

		newRnas = 0
		ntpsUsed = np.zeros(4)

		ntpsShape = self.bulkMoleculesPartition.ntps.counts().shape

		rnasCreated = np.zeros_like(self.bulkMoleculesPartition.rnas.counts())

		while enzLimit > 0:
			if not np.any(
					np.all(
						self.bulkMoleculesPartition.ntps.counts() > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not np.any(enzLimit > np.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if np.sum(self.rnaSynthProb[np.all(self.bulkMoleculesPartition.ntps.counts() > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if np.sum(self.rnaSynthProb[enzLimit > np.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdx = np.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]

			if np.any(self.bulkMoleculesPartition.ntps.counts() < self.rnaNtCounts[newIdx, :]):
				break

			if enzLimit < np.sum(self.rnaNtCounts[newIdx, :]):
				break

			enzLimit -= np.sum(self.rnaNtCounts[newIdx, :])

			# ntpsUsed += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)
			# self.metabolite.parentState.tcNtpUsage += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)

			self.bulkMoleculesPartition.ntps.countsDec(
				self.rnaNtCounts[newIdx, :].reshape(ntpsShape)
				)

			self.bulkMoleculesPartition.h2oMol.countDec(1)
			self.bulkMoleculesPartition.ppiMol.countInc(self.rnaLens[newIdx])
			self.bulkMoleculesPartition.hMol.countInc(1)

			rnasCreated[newIdx] += 1

			# Increment RNA
			# self.rna.counts[newIdx] += 1
			# newRnas += 1

		self.bulkMoleculesPartition.rnas.countsInc(rnasCreated)

		# print "%d" % enzLimit

#		print "Transcription newRnas: %d" % newRnas
#		print "Transcription ntpsUsed: %s" % str(ntpsUsed)
#		print "Transcription numActiveRnaps (total): %d (%d)" % (int(np.sum(ntpsUsed) / self.elngRate / self.timeStepSec), int(calcRnaps(self.enzyme.counts)))


def calcRnaps(countRpoA, countRpoB, countRpoC, countRpoD):
	return np.min([
		np.floor(countRpoA/2),
		countRpoB,
		countRpoC,
		countRpoD,
		])

_metIds = [
	"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
	"PPI[c]", "H2O[c]", "H[c]",
	]
