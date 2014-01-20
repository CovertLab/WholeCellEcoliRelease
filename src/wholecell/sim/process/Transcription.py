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
		
		# References to states
		self.metabolite = None
		self.rna = None
		self.enzyme = None

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

		mc = sim.getState("MoleculeCounts")

		# Metabolites
		self.metabolite = mc.addPartition(self, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], self.calcReqMetabolites)

		self.metaboliteView = mc.countsBulkViewNew([
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			])

		# self.metabolite.idx["ntps"] = self.metabolite.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[0]
		# self.metabolite.idx["ppi"] = self.metabolite.getIndex(["PPI[c]"])[0]
		# self.metabolite.idx["h2o"] = self.metabolite.getIndex(["H2O[c]"])[0]
		# self.metabolite.idx["h"] = self.metabolite.getIndex(["H[c]"])[0]

		self.metabolite.ntpView = self.metabolite.countsBulkViewNew(["ATP", "CTP", "GTP", "UTP"])
		self.metabolite.ppiMol = self.metabolite.molecule('PPI:mature', 'merged')
		self.metabolite.h2oMol = self.metabolite.molecule('H2O:mature', 'merged')
		self.metabolite.hMol = self.metabolite.molecule('H:mature', 'merged')

		# RNA
		self.rna = mc.addPartition(self, [x["id"] + ":nascent[c]" for x in kb.rnas], self.calcReqRna)
		self.rnaNtCounts = numpy.array([x["ntCount"] for x in kb.rnas])
		self.rnaLens = numpy.sum(self.rnaNtCounts, axis = 1)
		# self.rnaSynthProb = mc.rnaExp * (numpy.log(2) / self.cellCycleLength + 1 / numpy.array([x["halfLife"] for x in kb.rnas]))
		# self.rnaSynthProb /= numpy.sum(self.rnaSynthProb)

		# Enzymes
		# self.enzyme = sim.getState("MoleculeCounts").addPartition(self, ["RNAP70-CPLX:mature[c]"], self.calcReqEnzyme)
		# self.enzyme.idx["rnaPol"] = self.enzyme.getIndex(["RNAP70-CPLX:mature[c]"])[0]
		self.enzyme = mc.addPartition(self, [
			"EG10893-MONOMER", "RPOB-MONOMER", "RPOC-MONOMER", "RPOD-MONOMER"
			], self.calcReqEnzyme)
		# self.enzyme.idx["rpoA"] = self.enzyme.getIndex(["EG10893-MONOMER"])[0]
		# self.enzyme.idx["rpoB"] = self.enzyme.getIndex(["RPOB-MONOMER"])[0]
		# self.enzyme.idx["rpoC"] = self.enzyme.getIndex(["RPOC-MONOMER"])[0]
		# self.enzyme.idx["rpoD"] = self.enzyme.getIndex(["RPOD-MONOMER"])[0]

		self.enzyme.rpoAMol = self.enzyme.molecule('EG10893-MONOMER:mature', 'merged')
		self.enzyme.rpoBMol = self.enzyme.molecule('RPOB-MONOMER:mature', 'merged')
		self.enzyme.rpoCMol = self.enzyme.molecule('RPOC-MONOMER:mature', 'merged')
		self.enzyme.rpoDMol = self.enzyme.molecule('RPOD-MONOMER:mature', 'merged')

		self.rpoAMol = mc.molecule('EG10893-MONOMER:mature', 'c')
		self.rpoBMol = mc.molecule('RPOB-MONOMER:mature', 'c')
		self.rpoCMol = mc.molecule('RPOC-MONOMER:mature', 'c')
		self.rpoDMol = mc.molecule('RPOD-MONOMER:mature', 'c')


	def calcRnaps(self, countRpoA, countRpoB, countRpoC, countRpoD):
		# return numpy.min((	numpy.floor(numpy.sum(counts[self.enzyme.idx["rpoA"]]) / 2),
		# 					numpy.sum(counts[self.enzyme.idx["rpoB"]]),
		# 					numpy.sum(counts[self.enzyme.idx["rpoC"]]),
		# 					numpy.sum(counts[self.enzyme.idx["rpoD"]])
		# 				))

		return numpy.min([
			numpy.floor(countRpoA/2),
			countRpoB,
			countRpoC,
			countRpoD,
			])

	# Calculate needed metabolites
	def calcReqMetabolites(self, request):
		request.ntpView.countsBulkIs(
			numpy.min([
				self.calcRnaps(
					self.rpoAMol.countBulk(), self.rpoBMol.countBulk(),
					self.rpoCMol.countBulk(), self.rpoDMol.countBulk()
					) * self.elngRate * self.timeStepSec,
				4 * numpy.min(self.metaboliteView.countsBulk())
				])
			)


	# Calculate needed RNA
	def calcReqRna(self, request):
		request.countsBulkIs(0)


	# Calculate needed enzymes
	def calcReqEnzyme(self, request):
		request.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		enzLimit = numpy.min([
			self.calcRnaps(
				self.enzyme.rpoAMol.countBulk(), self.enzyme.rpoBMol.countBulk(),
				self.enzyme.rpoCMol.countBulk(), self.enzyme.rpoDMol.countBulk()
				) * self.elngRate * self.timeStepSec,
			1.1 * 4 * numpy.min(self.metabolite.countsBulk())
			])

		newRnas = 0
		ntpsUsed = numpy.zeros(4)

		while enzLimit > 0:
			print 'not implemented yet!'
			
			import ipdb
			ipdb.set_trace()

			if not numpy.any(
					numpy.all(
						self.metabolite.counts[self.metabolite.idx["ntps"]] > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not numpy.any(enzLimit > numpy.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if numpy.sum(self.rnaSynthProb[numpy.all(self.metabolite.counts[self.metabolite.idx["ntps"]] > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if numpy.sum(self.rnaSynthProb[enzLimit > numpy.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdx = numpy.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]

			if numpy.any(self.metabolite.counts[self.metabolite.idx["ntps"]] < self.rnaNtCounts[newIdx, :]):
				break

			if enzLimit < numpy.sum(self.rnaNtCounts[newIdx, :]):
				break

			enzLimit -= numpy.sum(self.rnaNtCounts[newIdx, :])

			self.metabolite.counts[self.metabolite.idx["ntps"]] -= self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)
			self.metabolite.counts[self.metabolite.idx["h2o"]] -= 1
			self.metabolite.counts[self.metabolite.idx["ppi"]] += self.rnaLens[newIdx]
			self.metabolite.counts[self.metabolite.idx["h"]] += 1
			ntpsUsed += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)
			self.metabolite.parentState.tcNtpUsage += self.rnaNtCounts[newIdx, :].reshape(self.metabolite.idx["ntps"].shape)

			# Increment RNA
			self.rna.counts[newIdx] += 1
			newRnas += 1
		# print "%d" % enzLimit

#		print "Transcription newRnas: %d" % newRnas
#		print "Transcription ntpsUsed: %s" % str(ntpsUsed)
#		print "Transcription numActiveRnaps (total): %d (%d)" % (int(numpy.sum(ntpsUsed) / self.elngRate / self.timeStepSec), int(self.calcRnaps(self.enzyme.counts)))