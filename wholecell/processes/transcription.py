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

		self.bulkCountsPartition.initialize(_metIds + rnaIds + enzIds)

		# Metabolites
		self.bulkCountsPartition.metabolites = self.bulkCountsPartition.countsBulkViewNew(_metIds)

		self.ntpView = mc.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.bulkCountsPartition.ntps = self.bulkCountsPartition.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.bulkCountsPartition.ppiMol = self.bulkCountsPartition.molecule('PPI[c]')
		self.bulkCountsPartition.h2oMol = self.bulkCountsPartition.molecule('H2O[c]')
		self.bulkCountsPartition.hMol = self.bulkCountsPartition.molecule('H[c]')

		# RNA
		self.bulkCountsPartition.rnas = self.bulkCountsPartition.countsBulkViewNew(rnaIds)
		self.rnaNtCounts = np.array([x["ntCount"] for x in kb.rnas])
		self.rnaLens = np.sum(self.rnaNtCounts, axis = 1)
		
		halflives = np.array([x["halfLife"] for x in kb.rnas])
		self.rnaSynthProb = mc.rnaExp * (np.log(2) / self.cellCycleLength + 1 / halflives)
		self.rnaSynthProb /= np.sum(self.rnaSynthProb)

		# Enzymes
		# self.enzyme = sim.states["BulkCounts"].addPartition(self, ["RNAP70-CPLX[c]"], self.calcReqEnzyme)
		# self.enzyme.idx["rnaPol"] = self.enzyme.getIndex(["RNAP70-CPLX[c]"])[0]
		self.bulkCountsPartition.enzymes = self.bulkCountsPartition.countsBulkViewNew(enzIds)
		self.bulkCountsPartition.rpoAMol = self.bulkCountsPartition.molecule('EG10893-MONOMER[c]')
		self.bulkCountsPartition.rpoBMol = self.bulkCountsPartition.molecule('RPOB-MONOMER[c]')
		self.bulkCountsPartition.rpoCMol = self.bulkCountsPartition.molecule('RPOC-MONOMER[c]')
		self.bulkCountsPartition.rpoDMol = self.bulkCountsPartition.molecule('RPOD-MONOMER[c]')

		self.rpoAMol = mc.molecule('EG10893-MONOMER[c]')
		self.rpoBMol = mc.molecule('RPOB-MONOMER[c]')
		self.rpoCMol = mc.molecule('RPOC-MONOMER[c]')
		self.rpoDMol = mc.molecule('RPOD-MONOMER[c]')


	def requestBulkCounts(self):
		self.bulkCountsPartition.ntps.countsBulkIs(
			np.min([
				calcRnaps(
					self.rpoAMol.countBulk(), self.rpoBMol.countBulk(),
					self.rpoCMol.countBulk(), self.rpoDMol.countBulk()
					) * self.elngRate * self.timeStepSec,
				4 * np.min(self.ntpView.countsBulk())
				])/4
			)

		self.bulkCountsPartition.h2oMol.countBulkIs(1)

		self.bulkCountsPartition.rnas.countsBulkIs(0)

		self.bulkCountsPartition.enzymes.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		enzLimit = np.min([
			calcRnaps(
				self.bulkCountsPartition.rpoAMol.countBulk(),
				self.bulkCountsPartition.rpoBMol.countBulk(),
				self.bulkCountsPartition.rpoCMol.countBulk(),
				self.bulkCountsPartition.rpoDMol.countBulk()
				) * self.elngRate * self.timeStepSec,
			1.1 * 4 * np.min(self.bulkCountsPartition.ntps.countsBulk())
			])

		newRnas = 0
		ntpsUsed = np.zeros(4)

		ntpsShape = self.bulkCountsPartition.ntps.countsBulk().shape

		rnasCreated = np.zeros_like(self.bulkCountsPartition.rnas.countsBulk())

		while enzLimit > 0:
			if not np.any(
					np.all(
						self.bulkCountsPartition.ntps.countsBulk() > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not np.any(enzLimit > np.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if np.sum(self.rnaSynthProb[np.all(self.bulkCountsPartition.ntps.countsBulk() > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if np.sum(self.rnaSynthProb[enzLimit > np.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdx = np.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]

			if np.any(self.bulkCountsPartition.ntps.countsBulk() < self.rnaNtCounts[newIdx, :]):
				break

			if enzLimit < np.sum(self.rnaNtCounts[newIdx, :]):
				break

			enzLimit -= np.sum(self.rnaNtCounts[newIdx, :])

			# ntpsUsed += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)
			# self.metabolite.parentState.tcNtpUsage += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)

			self.bulkCountsPartition.ntps.countsBulkDec(
				self.rnaNtCounts[newIdx, :].reshape(ntpsShape)
				)

			self.bulkCountsPartition.h2oMol.countBulkDec(1)
			self.bulkCountsPartition.ppiMol.countBulkInc(self.rnaLens[newIdx])
			self.bulkCountsPartition.hMol.countBulkInc(1)

			rnasCreated[newIdx] += 1

			# Increment RNA
			# self.rna.counts[newIdx] += 1
			# newRnas += 1

		self.bulkCountsPartition.rnas.countsBulkInc(rnasCreated)

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
