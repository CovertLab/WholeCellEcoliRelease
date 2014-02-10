#!/usr/bin/env python

"""
Transcription

Transcription sub-model. Encodes molecular simulation of macromolecular bacterial transcription

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.sim.process.Process

class Transcription(wholecell.sim.process.Process.Process):
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
		self.cellCycleLength = 1 * 3600		# s
		self.elngRate = 50					# nt/s
		self.rnaLens = None					# RNA lengths
		self.rnaNtCounts = None				# RNA nucleotide counts [nt x RNA] <-- TODO: Check this
		self.rnaSynthProb = None			# Relative RNA synthesis rates

		super(Transcription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Transcription, self).initialize(sim, kb)

		mc = sim.states["MoleculeCounts"]

		rnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas]

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.mcPartition = mc.setPartition(self, _metIds + rnaIds + enzIds)

		# Metabolites
		self.mcPartition.metabolites = self.mcPartition.countsBulkViewNew(_metIds)

		self.ntpView = mc.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.mcPartition.ntps = self.mcPartition.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])

		self.mcPartition.ppiMol = self.mcPartition.molecule('PPI[c]')
		self.mcPartition.h2oMol = self.mcPartition.molecule('H2O[c]')
		self.mcPartition.hMol = self.mcPartition.molecule('H[c]')

		# RNA
		self.mcPartition.rnas = self.mcPartition.countsBulkViewNew(rnaIds)
		self.rnaNtCounts = numpy.array([x["ntCount"] for x in kb.rnas])
		self.rnaLens = numpy.sum(self.rnaNtCounts, axis = 1)
		self.rnaSynthProb = mc.rnaExp * (numpy.log(2) / self.cellCycleLength + 1 / numpy.array([x["halfLife"] for x in kb.rnas]))
		self.rnaSynthProb /= numpy.sum(self.rnaSynthProb)

		# Enzymes
		# self.enzyme = sim.states["MoleculeCounts"].addPartition(self, ["RNAP70-CPLX[c]"], self.calcReqEnzyme)
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


	def calcRnaps(self, countRpoA, countRpoB, countRpoC, countRpoD):
		return numpy.min([
			numpy.floor(countRpoA/2),
			countRpoB,
			countRpoC,
			countRpoD,
			])

	def requestMoleculeCounts(self):
		self.mcPartition.ntps.countsBulkIs(
			numpy.min([
				self.calcRnaps(
					self.rpoAMol.countBulk(), self.rpoBMol.countBulk(),
					self.rpoCMol.countBulk(), self.rpoDMol.countBulk()
					) * self.elngRate * self.timeStepSec,
				4 * numpy.min(self.ntpView.countsBulk())
				])/4
			)

		self.mcPartition.h2oMol.countBulkIs(1)

		self.mcPartition.rnas.countsBulkIs(0)

		self.mcPartition.enzymes.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		enzLimit = numpy.min([
			self.calcRnaps(
				self.mcPartition.rpoAMol.countBulk(),
				self.mcPartition.rpoBMol.countBulk(),
				self.mcPartition.rpoCMol.countBulk(),
				self.mcPartition.rpoDMol.countBulk()
				) * self.elngRate * self.timeStepSec,
			1.1 * 4 * numpy.min(self.mcPartition.ntps.countsBulk())
			])

		newRnas = 0
		ntpsUsed = numpy.zeros(4)

		ntpsShape = self.mcPartition.ntps.countsBulk().shape

		rnasCreated = numpy.zeros_like(self.mcPartition.countsBulk())

		while enzLimit > 0:
			if not numpy.any(
					numpy.all(
						self.mcPartition.ntps.countsBulk() > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not numpy.any(enzLimit > numpy.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if numpy.sum(self.rnaSynthProb[numpy.all(self.mcPartition.ntps.countsBulk() > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if numpy.sum(self.rnaSynthProb[enzLimit > numpy.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdx = numpy.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]

			if numpy.any(self.mcPartition.ntps.countsBulk() < self.rnaNtCounts[newIdx, :]):
				break

			if enzLimit < numpy.sum(self.rnaNtCounts[newIdx, :]):
				break

			enzLimit -= numpy.sum(self.rnaNtCounts[newIdx, :])

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

		self.mcPartition.countsBulkInc(rnasCreated)

		# print "%d" % enzLimit

#		print "Transcription newRnas: %d" % newRnas
#		print "Transcription ntpsUsed: %s" % str(ntpsUsed)
#		print "Transcription numActiveRnaps (total): %d (%d)" % (int(numpy.sum(ntpsUsed) / self.elngRate / self.timeStepSec), int(self.calcRnaps(self.enzyme.counts)))

_metIds = [
	"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
	"PPI[c]", "H2O[c]", "H[c]",
	]
