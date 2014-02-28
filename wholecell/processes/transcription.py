#!/usr/bin/env python

"""
Transcription

Transcription sub-model. Encodes molecular simulation of macromolecular bacterial transcription

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

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

		mc = sim.states["BulkCounts"]

		rnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas]

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.mcPartition.initialize(_metIds + rnaIds + enzIds)

		# Metabolites
		self.mcPartition.metabolites = self.mcPartition.countsBulkViewNew(_metIds)

		self.ntpView = mc.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.mcPartition.ntps = self.mcPartition.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.mcPartition.ppiMol = self.mcPartition.molecule('PPI[c]')
		self.mcPartition.h2oMol = self.mcPartition.molecule('H2O[c]')
		self.mcPartition.hMol = self.mcPartition.molecule('H[c]')

		# RNA
		self.mcPartition.rnas = self.mcPartition.countsBulkViewNew(rnaIds)
		self.rnaNtCounts = np.array([x["ntCount"] for x in kb.rnas])
		self.rnaLens = np.sum(self.rnaNtCounts, axis = 1)
		
		halflives = np.array([x["halfLife"] for x in kb.rnas])
		self.rnaSynthProb = mc.rnaExp * (np.log(2) / self.cellCycleLength + 1 / halflives)
		self.rnaSynthProb /= np.sum(self.rnaSynthProb)

		# Enzymes
		# self.enzyme = sim.states["BulkCounts"].addPartition(self, ["RNAP70-CPLX[c]"], self.calcReqEnzyme)
		# self.enzyme.idx["rnaPol"] = self.enzyme.getIndex(["RNAP70-CPLX[c]"])[0]
		self.mcPartition.enzymes = self.mcPartition.countsBulkViewNew(enzIds)
		self.mcPartition.rpoAMol = self.mcPartition.molecule('EG10893-MONOMER[c]')
		self.mcPartition.rpoBMol = self.mcPartition.molecule('RPOB-MONOMER[c]')
		self.mcPartition.rpoCMol = self.mcPartition.molecule('RPOC-MONOMER[c]')
		self.mcPartition.rpoDMol = self.mcPartition.molecule('RPOD-MONOMER[c]')

		self.rpoAMol = mc.molecule('EG10893-MONOMER[c]')
		self.rpoBMol = mc.molecule('RPOB-MONOMER[c]')
		self.rpoCMol = mc.molecule('RPOC-MONOMER[c]')
		self.rpoDMol = mc.molecule('RPOD-MONOMER[c]')


	def requestBulkCounts(self):
		self.mcPartition.ntps.countsBulkIs(
			np.min([
				calcRnaps(
					self.rpoAMol.countBulk(), self.rpoBMol.countBulk(),
					self.rpoCMol.countBulk(), self.rpoDMol.countBulk()
					) * self.elngRate * self.timeStepSec,
				4 * np.min(self.ntpView.countsBulk())
				])/4
			)

		self.mcPartition.h2oMol.countBulkIs(1)

		self.mcPartition.rnas.countsBulkIs(0)

		self.mcPartition.enzymes.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		enzLimit = np.min([
			calcRnaps(
				self.mcPartition.rpoAMol.countBulk(),
				self.mcPartition.rpoBMol.countBulk(),
				self.mcPartition.rpoCMol.countBulk(),
				self.mcPartition.rpoDMol.countBulk()
				) * self.elngRate * self.timeStepSec,
			1.1 * 4 * np.min(self.mcPartition.ntps.countsBulk())
			])

		newRnas = 0
		ntpsUsed = np.zeros(4)

		ntpsShape = self.mcPartition.ntps.countsBulk().shape

		rnasCreated = np.zeros_like(self.mcPartition.rnas.countsBulk())

		while enzLimit > 0:
			if not np.any(
					np.all(
						self.mcPartition.ntps.countsBulk() > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not np.any(enzLimit > np.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if np.sum(self.rnaSynthProb[np.all(self.mcPartition.ntps.countsBulk() > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if np.sum(self.rnaSynthProb[enzLimit > np.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdx = np.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]

			if np.any(self.mcPartition.ntps.countsBulk() < self.rnaNtCounts[newIdx, :]):
				break

			if enzLimit < np.sum(self.rnaNtCounts[newIdx, :]):
				break

			enzLimit -= np.sum(self.rnaNtCounts[newIdx, :])

			# ntpsUsed += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)
			# self.metabolite.parentState.tcNtpUsage += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)

			self.mcPartition.ntps.countsBulkDec(
				self.rnaNtCounts[newIdx, :].reshape(ntpsShape)
				)

			self.mcPartition.h2oMol.countBulkDec(1)
			self.mcPartition.ppiMol.countBulkInc(self.rnaLens[newIdx])
			self.mcPartition.hMol.countBulkInc(1)

			rnasCreated[newIdx] += 1

			# Increment RNA
			# self.rna.counts[newIdx] += 1
			# newRnas += 1

		self.mcPartition.rnas.countsBulkInc(rnasCreated)

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
