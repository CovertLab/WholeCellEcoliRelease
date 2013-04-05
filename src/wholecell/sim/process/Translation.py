#!/usr/bin/env python

"""
Translation

Translation sub-model. Encodes molecular simulation of macromolecular bacterial translation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.sim.process.Process

class Translation(wholecell.sim.process.Process.Process):
	""" Translation """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "Translation",
			"name": "Translation"
		}
		
		# References to states
		self.metabolite = None
		self.mrna = None
		self.protein = None
		self.enzyme = None

		# Constants
		self.elngRate = 16			# AA/s
		self.proteinAaCounts = None	# Protein amino acid counts [AA x protein] <-- TODO: Check this
		self.proteinLens = None		# Protein lengths

		super(Translation, self).__init__()

	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		super(Translation, self).initialize(sim, kb)

		# Metabolites
		self.metabolite = sim.getState("MoleculeCounts").addPartition(self, [
				"ALA[c]", "ARG[c]", "ASN[c]", "ASP[c]", "CYS[c]", "GLU[c]", "GLN[c]", "GLY[c]", "HIS[c]", "ILE[c]",  "LEU[c]",
				"LYS[c]", "MET[c]", "PHE[c]", "PRO[c]", "SER[c]", "THR[c]", "TRP[c]", "TYR[c]", "VAL[c]",
				"FMET[c]",
				"GTP[c]", "GDP[c]", "PI[c]",  "H2O[c]", "H[c]"
			], self.calcReqMetabolites)

		self.metabolite.idx["aas"] = self.metabolite.getIndex([
			"ALA[c]", "ARG[c]", "ASN[c]", "ASP[c]", "CYS[c]", "GLU[c]", "GLN[c]", "GLY[c]", "HIS[c]", "ILE[c]",  "LEU[c]",
			"LYS[c]", "MET[c]", "PHE[c]", "PRO[c]", "SER[c]", "THR[c]", "TRP[c]", "TYR[c]", "VAL[c]",
			])[0]
		self.metabolite.idx["gtp"] = self.metabolite.getIndex(["GTP[c]"])[0]
		self.metabolite.idx["gdp"] = self.metabolite.getIndex(["GDP[c]"])[0]
		self.metabolite.idx["pi"] = self.metabolite.getIndex(["PI[c]"])[0]
		self.metabolite.idx["h2o"] = self.metabolite.getIndex(["H2O[c]"])[0]
		self.metabolite.idx["h"] = self.metabolite.getIndex(["H[c]"])[0]

		# mRNA, protein monomer
		mrnas = [x for x in kb.rnas if x["type"] == "mRNA"]
		monomers = [x for x in kb.proteins if x["monomer"] == True]
		self.mrna = sim.getState("MoleculeCounts").addPartition(self,[x["id"] + ":mature[c]" for x in mrnas], self.calcReqMrna)
		self.mrna = sim.getState("MoleculeCounts").addPartition(self,[x["monomerId"] + ":nascent[c]" for x in mrnas], self.calcReqProtein)
		self.proteinAaCounts = numpy.array([x["aaCount"] for x in monomers])
		self.proteinLens = numpy.sum(self.proteinAaCounts, axis = 1)

		# Enzymes
		self.enzyme = sim.getState("MoleculeCounts").addPartition(self, ["RIBOSOME_70S:mature[c]"], self.calcReqEnzyme)
		self.enzyme.idx["ribosome70S"] = self.enzyme.getIndex("RIBOSOME_70S:mature[c]")[0]

	# Calculate needed metabolites
	def calcReqMetabolites(self):
		val = numpy.zeros(self.metabolite.fullCounts.shape)

		elng = numpy.min([
			self.enzyme.fullCounts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec,		# Polymerization by all available ribosomes
			numpy.sum(self.metabolite.fullCounts[self.metabolite.idx["aas"]])								# Amino acid limitation
			])
		val[self.metabolite.idx["aas"]] = elng / self.metabolite.idx["aas"].size
		val[[self.metabolite.idx["gtp"], self.metabolite.idx["h2o"]]] = 4 * 2 * elng
		return val

	# Calculate needed mRNA
	def calcReqMrna(self):
		return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate needed protein monomers
	def calcReqProtein(self):
		return numpy.zeros(self.enzyme.fullCounts.shape)

	# Calculate needed enzymes
	def calcReqEnzyme(self):
		return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		if not numpy.any(self.mrna.counts):
			return

		# Total synthesis rate
		proteinSynthProb = self.mrna.counts / numpy.sum(self.mrna.counts)
		totRate = 1 / numpy.dot(self.proteinLens, self.proteinSynthProb) * numpy.min([					# Normalize by average protein length
			numpy.sum(self.metabolite.counts[self.metabolite.idx["aas"]]),								# Amino acid limitation
			numpy.sum(self.metabolite.counts[self.metabolite.idx["gtp"]]) / 2,							# GTP (energy) limitation
			self.enzyme.counts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec		# Ribosome capacity
			])

		# Gillespie-like algorithm
		t = 0
		while True:
			# Choose time step
			t -= numpy.log(self.randStream.rand()) / totRate
			if t > self.timeStepSec:
				break

			# Check if sufficient metabolic resources to make protein
			newIdx = numpy.where(self.randStream.mnrnd(1, self.proteinSynthProb))[0]
			if \
				numpy.any(self.proteinAaCounts[newIdx, :] > self.metabolite.counts[self.metabolite.idx["aas"]]) or \
				2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["gtp"]] or \
				2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["h2o"]]:
					break

			# Update metabolites
			self.metabolite.counts[self.metabolite.idx["aas"]] -= self.proteinAaCounts[newIdx, :]
			self.metabolite.counts[self.metabolite.idx["h2o"]] += self.proteinLens[newIdx] - 1

			self.metabolite.counts[self.metabolite.idx["gtp"]] -= 2 * self.proteinLens[newIdx]
			self.metabolite.counts[self.metabolite.idx["h2o"]] -= 2 * self.proteinLens[newIdx]
			self.metabolite.counts[self.metabolite.idx["gdp"]] += 2 * self.proteinLens[newIdx]
			self.metabolite.counts[self.metabolite.idx["pi"]] += 2 * self.proteinLens[newIdx]
			self.metabolite.counts[self.metabolite.idx["h"]] += 2 * self.proteinLens[newIdx]

			# Increment protein monomer
			self.protein.counts[newIdx] += 1