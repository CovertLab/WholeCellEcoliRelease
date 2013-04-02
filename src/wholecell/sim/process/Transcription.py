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

	meta = {
		"id": "Transcription",
		"name": "Transcription"
	}

	# Constructor
	def __init__(self):
		# References to states
		self.metabolite = None
		self.rna = None
		self.enzyme = None

		# Constants
		self.cellCycleLength = 9 * 3600		# s
		self.elngRate = 50					# nt/s
		self.rnaLens = None					# RNA lengths
		self.rnaNtCounts = None				# RNA nucleotide counts [nt x RNA] <-- TODO: Check this
		self.rnaSynthProb = None			# Relative RNA synthesis rates

		super(Transcription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Transcription, self).initialize(sim, kb)

		# Metabolites
		self.metabolite = sim.getState("MoleculeCounts").addPartition(self, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], self.calcReqMetabolites)

		self.metabolite.idx["ntps"] = self.metabolite.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.metabolite.idx["ppi"] = self.metabolite.getIndex(["PPI[c]"])
		self.metabolite.idx["h2o"] = self.metabolite.getIndex(["H2O[c]"])
		self.metabolite.idx["h"] = self.metabolite.getIndex(["H[c]"])

		# RNA
		self.rna = sim.getState("MoleculeCounts").addPartition(self, [x["id"] + ":nascent[c]" for x in kb.rnas], self.calcReqRna)
		self.rnaNtCounts = numpy.array([x["ntCount"] for x in kb.rnas])
		self.rnaLens = numpy.sum(self.rnaNtCounts, axis = 1)
		self.rnaSynthProb =  numpy.array([x["exp"] for x in kb.rnas]) * (numpy.log(2) / self.cellCycleLength + 1 / numpy.array([x["halfLife"] for x in kb.rnas]))
		self.rnaSynthProb /= numpy.sum(self.rnaSynthProb)

		# Enzymes
		self.enzyme = sim.getState("MoleculeCounts").addPartition(self, ["RNA_POLYMERASE:mature[c]"], self.calcReqEnzyme)
		self.enzyme.idx["rnaPol"] = self.enzyme.getIndex(["RNA_POLYMERASE:mature[c]"])

	# Calculate needed metabolites
	def calcReqMetabolites(self):
		val = numpy.zeros(self.metabolite.fullCounts.shape)

		val[self.metabolite.idx["ntps"]] = numpy.min([
			self.enzyme.fullCounts[self.enzyme.idx["rnaPol"]] * self.elngRate * self.timeStepSec,		# Polymerization by all available RNA Polymerases
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
		# Total synthesis rate
		totRate = numpy.minimum(
			numpy.sum(self.metabolite.counts[self.metabolite.idx["ntps"]]),								# NTP Limitation
			self.enzyme.counts[self.enzyme.idx["rnaPol"]] * self.elngRate * self.timeStepSec			# Polymerization by all available RNA Polymerases
			) / numpy.dot(self.rnaLens, self.rnaSynthProb)												# Normalize by expected NTP usage

		# Gillespie-like algorithm
		t = 0
		while True:
			# Choose time step
			t -= numpy.log(self.randStream.rand()) / totRate
			if t > self.timeStepSec:
				break

			# Check if sufficient metabolic resources to make RNA
			newIdx = numpy.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]
			if \
				numpy.any(self.rnaNtCounts[newIdx, :] > self.metabolite.counts[self.metabolite.idx["ntps"]]) or \
				self.metabolite.counts[self.metabolite.idx["h2o"]] < 1:
					break

			# Update metabolites
			self.metabolite.counts[self.metabolite.idx["ntps"]] -= self.rnaNtCounts[newIdx, :]
			self.metabolite.counts[self.metabolite.idx["h2o"]] -= 1
			self.metabolite.counts[self.metabolite.idx["ppi"]] += self.rnaLens(newIdx)
			self.metabolite.counts[self.metabolite.idx["h"]] += 1

			# Increment RNA
			self.rna.counts[newIdx] += 1