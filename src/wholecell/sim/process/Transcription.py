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

		self.metabolite.idx["ntps"] = self.metabolite.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[0]
		self.metabolite.idx["ppi"] = self.metabolite.getIndex(["PPI[c]"])[0]
		self.metabolite.idx["h2o"] = self.metabolite.getIndex(["H2O[c]"])[0]
		self.metabolite.idx["h"] = self.metabolite.getIndex(["H[c]"])[0]

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
		self.enzyme.idx["rpoA"] = self.enzyme.getIndex(["EG10893-MONOMER"])[0]
		self.enzyme.idx["rpoB"] = self.enzyme.getIndex(["RPOB-MONOMER"])[0]
		self.enzyme.idx["rpoC"] = self.enzyme.getIndex(["RPOC-MONOMER"])[0]
		self.enzyme.idx["rpoD"] = self.enzyme.getIndex(["RPOD-MONOMER"])[0]

	def calcRnaps(self, counts):
		return numpy.min((	numpy.floor(numpy.sum(counts[self.enzyme.idx["rpoA"]]) / 2),
							numpy.sum(counts[self.enzyme.idx["rpoB"]]),
							numpy.sum(counts[self.enzyme.idx["rpoC"]]),
							numpy.sum(counts[self.enzyme.idx["rpoD"]])
						))

	# Calculate needed metabolites
	def calcReqMetabolites(self):
		val = numpy.zeros(self.metabolite.fullCounts.shape)

		val[self.metabolite.idx["ntps"]] = numpy.min([
			self.calcRnaps(self.enzyme.fullCounts) * self.elngRate * self.timeStepSec,
			# self.enzyme.fullCounts[self.enzyme.idx["rnaPol"]] * self.elngRate * self.timeStepSec,		# Polymerization by all available RNA Polymerases
			4 * numpy.min(self.metabolite.fullCounts[self.metabolite.idx["ntps"]])						# Limited by scarcest NTP
			]) / 4

		val[self.metabolite.idx["h2o"]] = 1
		return val

	# Calculate needed RNA
	def calcReqRna(self):
		return numpy.zeros(self.rna.fullCounts.shape)

	# Calculate needed enzymes
	def calcReqEnzyme(self):
		return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		print "TRANSCRIPTION"
		enzLimit = numpy.min([
			self.calcRnaps(self.enzyme.counts) * self.elngRate * self.timeStepSec,
			1.1 * 4 * numpy.min(self.metabolite.counts[self.metabolite.idx["ntps"]])
			])

		print "Transcription enzLimit: %0.3f" % (enzLimit)
		print "Transcription ntps: %s" % str(self.metabolite.counts[self.metabolite.idx["ntps"]])
		import sys
		sys.stdout.flush()

		newRnas = 0
		ntpsUsed = numpy.zeros(4)
		self.metabolite.parentState.tcNtpUsage = numpy.zeros(4)

		while enzLimit > 0:
			if not numpy.any(numpy.all(self.metabolite.counts[self.metabolite.idx["ntps"]] > self.rnaNtCounts, axis = 1)):
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
				continue

			if enzLimit < numpy.sum(self.rnaNtCounts[newIdx, :]):
				continue

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

		print "Transcription newRnas: %d" % newRnas
		print "Transcription ntpsUsed: %s" % str(ntpsUsed)
		print "Transcription numActiveRnaps: %d" % int(numpy.sum(ntpsUsed) / self.elngRate / self.timeStepSec)