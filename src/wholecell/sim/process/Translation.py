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
				"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
				"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
				#"FMET[c]", # TODO: Re-add
				# "GTP[c]", "GDP[c]", "PI[c]",  "H2O[c]", "H[c]"
				"ATP[c]", "ADP[c]", "PI[c]",  "H2O[c]", "H[c]"
			], self.calcReqMetabolites)

		self.metabolite.idx["aas"] = self.metabolite.getIndex([
			"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
			"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
			])[0]
		self.metabolite.idx["atp"] = self.metabolite.getIndex(["ATP[c]"])[0]
		self.metabolite.idx["adp"] = self.metabolite.getIndex(["ADP[c]"])[0]
		self.metabolite.idx["pi"] = self.metabolite.getIndex(["PI[c]"])[0]
		self.metabolite.idx["h2o"] = self.metabolite.getIndex(["H2O[c]"])[0]
		self.metabolite.idx["h"] = self.metabolite.getIndex(["H[c]"])[0]

		# mRNA, protein monomer
		mrnas = [x for x in kb.rnas if x["monomerId"] != None]
		monomers = [x for x in kb.proteins if x["monomer"] == True and x["modifiedForm"] == False]
		self.mrna = sim.getState("MoleculeCounts").addPartition(self,[x["id"] + ":mature[c]" for x in mrnas], self.calcReqMrna)
		self.protein = sim.getState("MoleculeCounts").addPartition(self,[x["monomerId"] + ":nascent[c]" for x in mrnas], self.calcReqProtein)
		self.proteinAaCounts = numpy.array([x["aaCount"] for x in monomers])
		self.proteinLens = numpy.sum(self.proteinAaCounts, axis = 1)

		# Enzymes
		# TODO: We really want all the associated riboproteins as well (ie we want the complexes)
		self.enzyme = sim.getState("MoleculeCounts").addPartition(self, [
			"RRLA-RRNA:mature[c]", "RRLB-RRNA:mature[c]", "RRLC-RRNA:mature[c]", "RRLD-RRNA:mature[c]", "RRLE-RRNA:mature[c]", "RRLG-RRNA:mature[c]", "RRLH-RRNA:mature[c]",
			"RRSA-RRNA:mature[c]", "RRSB-RRNA:mature[c]", "RRSC-RRNA:mature[c]", "RRSD-RRNA:mature[c]", "RRSE-RRNA:mature[c]", "RRSG-RRNA:mature[c]", "RRSH-RRNA:mature[c]",
			"RRFA-RRNA:mature[c]", "RRFB-RRNA:mature[c]", "RRFC-RRNA:mature[c]", "RRFD-RRNA:mature[c]", "RRFE-RRNA:mature[c]", "RRFF-RRNA:mature[c]", "RRFG-RRNA:mature[c]", "RRFH-RRNA:mature[c]"
			], self.calcReqEnzyme)
		self.enzyme.idx["23S"] = self.enzyme.getIndex([
			"RRLA-RRNA:mature[c]", "RRLB-RRNA:mature[c]", "RRLC-RRNA:mature[c]", "RRLD-RRNA:mature[c]", "RRLE-RRNA:mature[c]", "RRLG-RRNA:mature[c]", "RRLH-RRNA:mature[c]"
			])[0]
		self.enzyme.idx["16S"] = self.enzyme.getIndex([
			"RRSA-RRNA:mature[c]", "RRSB-RRNA:mature[c]", "RRSC-RRNA:mature[c]", "RRSD-RRNA:mature[c]", "RRSE-RRNA:mature[c]", "RRSG-RRNA:mature[c]", "RRSH-RRNA:mature[c]",
			])[0]
		self.enzyme.idx["5S"] = self.enzyme.getIndex([
			"RRFB-RRNA:mature[c]", "RRFC-RRNA:mature[c]", "RRFD-RRNA:mature[c]", "RRFE-RRNA:mature[c]", "RRFF-RRNA:mature[c]", "RRFG-RRNA:mature[c]", "RRFH-RRNA:mature[c]"
			])[0]

	def calcRibosomes(self, counts):
		return numpy.min((numpy.sum(counts[self.enzyme.idx["23S"]]), numpy.sum(counts[self.enzyme.idx["16S"]]), numpy.sum(counts[self.enzyme.idx["5S"]])))

	# Calculate needed metabolites
	def calcReqMetabolites(self):
		val = numpy.zeros(self.metabolite.fullCounts.shape)

		elng = numpy.min([
			self.calcRibosomes(self.enzyme.fullCounts) * self.elngRate * self.timeStepSec,
			# self.enzyme.fullCounts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec,		# Polymerization by all available ribosomes
			numpy.sum(self.metabolite.fullCounts[self.metabolite.idx["aas"]])								# Amino acid limitation
			])
		val[self.metabolite.idx["aas"]] = elng / self.metabolite.idx["aas"].size
		# val[numpy.array([self.metabolite.idx["atp"], self.metabolite.idx["h2o"]])] = 4 * 2 * elng
		return val

	# Calculate needed mRNA
	def calcReqMrna(self):
		return numpy.ones(self.mrna.fullCounts.shape)

	# Calculate needed protein monomers
	def calcReqProtein(self):
		return numpy.zeros(self.protein.fullCounts.shape)

	# Calculate needed enzymes
	def calcReqEnzyme(self):
		return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		if not numpy.any(self.mrna.counts):
			return

		# Total synthesis rate
		proteinSynthProb = self.mrna.counts / numpy.sum(self.mrna.counts)
		totRate = 1 / numpy.dot(self.proteinLens, proteinSynthProb) * numpy.min([					# Normalize by average protein length
			numpy.sum(self.metabolite.counts[self.metabolite.idx["aas"]]),								# Amino acid limitation
			# numpy.sum(self.metabolite.counts[self.metabolite.idx["atp"]]) / 2,							# GTP (energy) limitation
			self.calcRibosomes(self.enzyme.counts) * self.elngRate * self.timeStepSec
			# self.enzyme.counts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec		# Ribosome capacity
			])

		print "Translation totRate: %0.3f" % (totRate)
		newProts = 0

		# Gillespie-like algorithm
		t = 0
		while True:
			# Choose time step
			t -= numpy.log(self.randStream.rand()) / totRate
			if t > self.timeStepSec:
				break

			# Check if sufficient metabolic resources to make protein
			newIdx = numpy.where(self.randStream.mnrnd(1, proteinSynthProb))[0]
			if \
				numpy.any(self.proteinAaCounts[newIdx, :] > self.metabolite.counts[self.metabolite.idx["aas"]]):
				# numpy.any(self.proteinAaCounts[newIdx, :] > self.metabolite.counts[self.metabolite.idx["aas"]]) or \
				# 2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["atp"]] or \
				# 2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["h2o"]]:
					break

			# Update metabolites
			self.metabolite.counts[self.metabolite.idx["aas"]] -= self.proteinAaCounts[newIdx, :].reshape(-1)
			# self.metabolite.counts[self.metabolite.idx["h2o"]] += self.proteinLens[newIdx] - 1

			# self.metabolite.counts[self.metabolite.idx["atp"]] -= 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["h2o"]] -= 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["adp"]] += 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["pi"]] += 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["h"]] += 2 * self.proteinLens[newIdx]

			# Increment protein monomer
			self.protein.counts[newIdx] += 1
			newProts += 1
		print "Translation newProts: %d" % newProts