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

		mc = sim.getState("MoleculeCounts")

		# Metabolites
		self.metabolite = mc.addPartition(self, _metIDs, self.calcReqMetabolites)

		self.metabolite.aas = self.metabolite.countsBulkViewNew(_aaIDs)

		self.n_aas = len(_aaIDs)

		self.metabolite.aasNotSec = self.metabolite.countsBulkViewNew(_aaNotSecIDs)

		self.aas = mc.countsBulkViewNew(_aaIDs)

		self.metabolite.atpMol = self.metabolite.molecule('ATP:mature', 'merged')
		self.metabolite.adpMol = self.metabolite.molecule('ADP:mature', 'merged')
		self.metabolite.piMol = self.metabolite.molecule('PI:mature', 'merged')
		self.metabolite.h2oMol = self.metabolite.molecule('H2O:mature', 'merged')
		self.metabolite.hMol = self.metabolite.molecule('H:mature', 'merged')

		# mRNA, protein monomer
		mrnas = [x for x in kb.rnas if x["monomerId"] != None]
		monomers = [x for x in kb.proteins if len(x["composition"]) == 0 and x["unmodifiedForm"] == None]
		self.mrna = mc.addPartition(self,[x["id"] + ":mature[c]" for x in mrnas], self.calcReqMrna)
		self.protein = mc.addPartition(self,[x["monomerId"] + ":nascent[c]" for x in mrnas], self.calcReqProtein)
		self.proteinAaCounts = numpy.array([x["aaCount"] for x in monomers])
		self.proteinLens = numpy.sum(self.proteinAaCounts, axis = 1)

		# Enzymes
		# TODO: We really want all the associated riboproteins as well (ie we want the complexes)
		self.enzyme = mc.addPartition(self, _enzIDs, self.calcReqEnzyme)

		self.enzyme.ribosome23S = self.enzyme.countsBulkViewNew(_rib23S_IDs)

		self.enzyme.ribosome16S = self.enzyme.countsBulkViewNew(_rib16S_IDs)

		self.enzyme.ribosome5S = self.enzyme.countsBulkViewNew(_rib5S_IDs)

		self.ribosome23S = mc.countsBulkViewNew(_rib23S_IDs)

		self.ribosome16S = mc.countsBulkViewNew(_rib16S_IDs)

		self.ribosome5S = mc.countsBulkViewNew(_rib5S_IDs)


	def calcRibosomes(self, counts23S, counts16S, counts5S):
		return numpy.min([
			numpy.sum(counts23S),
			numpy.sum(counts16S),
			numpy.sum(counts5S)
			])


	# Calculate needed metabolites
	def calcReqMetabolites(self, request):
		# val = numpy.zeros(self.metabolite.fullCounts.shape)

		# elng = numpy.min([
		# 	self.calcRibosomes(self.enzyme.fullCounts) * self.elngRate * self.timeStepSec,
		# 	# self.enzyme.fullCounts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec,		# Polymerization by all available ribosomes
		# 	numpy.sum(self.metabolite.fullCounts[self.metabolite.idx["aas"]])								# Amino acid limitation
		# 	])
		# val[self.metabolite.idx["aas"]] = elng / self.metabolite.idx["aas"].size
		# # val[numpy.array([self.metabolite.idx["atp"], self.metabolite.idx["h2o"]])] = 4 * 2 * elng
		# return val

		ribs = self.calcRibosomes(
			self.ribosome23S.countsBulk(),
			self.ribosome16S.countsBulk(),
			self.ribosome5S.countsBulk()
			)

		elng = numpy.min([
			ribs * self.elngRate * self.timeStepSec,
			numpy.sum(self.aas.countsBulk())
			])

		request.aas.countsBulkIs(elng / self.n_aas)


	# Calculate needed mRNA
	def calcReqMrna(self, request):
		request.countsBulkIs(1)


	# Calculate needed protein monomers
	def calcReqProtein(self, request):
		request.countsBulkIs(0)


	# Calculate needed enzymes
	def calcReqEnzyme(self, request):
		request.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		if not numpy.any(self.mrna.countsBulk()):
			return

		print 'not implemented yet!'
		
		import ipdb
		ipdb.set_trace()

		# print "nRibosomes: %d" % int(self.calcRibosomes(self.enzyme.counts))
		# Total synthesis rate
		proteinSynthProb = self.mrna.counts / numpy.sum(self.mrna.counts)
		totRate = 1 / numpy.dot(self.proteinLens, proteinSynthProb) * numpy.min([					# Normalize by average protein length
			numpy.sum(self.metabolite.counts[self.metabolite.idx["aas"]]),								# Amino acid limitation
			# numpy.sum(self.metabolite.counts[self.metabolite.idx["atp"]]) / 2,							# GTP (energy) limitation
			self.calcRibosomes(self.enzyme.counts) * self.elngRate * self.timeStepSec
			# self.enzyme.counts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec		# Ribosome capacity
			])

		# print "Translation totRate: %0.3f" % (totRate)
		newProts = 0
		aasUsed = numpy.zeros(21)

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
			aasUsed += self.proteinAaCounts[newIdx, :].reshape(-1)
		self.aasUsed = aasUsed
		# print "Translation newProts: %d" % newProts
		# print "Translation aasUsed: %s" % str(aasUsed)
		# print "Translation numActiveRibs (total): %d (%d)" % (int(numpy.sum(aasUsed) / self.elngRate / self.timeStepSec), int(self.calcRibosomes(self.enzyme.counts)))

_metIDs = ["ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]",
	"GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
	"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]",
	"THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	#"FMET[c]", # TODO: Re-add
	#"GTP[c]", "GDP[c]", "PI[c]",  "H2O[c]", "H[c]"
	"ATP[c]", "ADP[c]", "PI[c]",  "H2O[c]", "H[c]"]

_aaIDs = [
	"ALA-L", "ARG-L", "ASN-L", "ASP-L", "CYS-L", "GLU-L", "GLN-L", "GLY", "HIS-L", "ILE-L",  "LEU-L",
	"LYS-L", "MET-L", "PHE-L", "PRO-L", "SELNP", "SER-L", "THR-L", "TRP-L", "TYR-L", "VAL-L",
	]

_aaNotSecIDs = [
	"ALA-L", "ARG-L", "ASN-L", "ASP-L", "CYS-L", "GLU-L", "GLN-L", "GLY", "HIS-L", "ILE-L",
	"LEU-L", "LYS-L", "MET-L", "PHE-L", "PRO-L", "SER-L", "THR-L", "TRP-L", "TYR-L", "VAL-L",
	]

_enzIDs = [
	"RRLA-RRNA:mature[c]", "RRLB-RRNA:mature[c]", "RRLC-RRNA:mature[c]", "RRLD-RRNA:mature[c]", "RRLE-RRNA:mature[c]", "RRLG-RRNA:mature[c]", "RRLH-RRNA:mature[c]",
	"RRSA-RRNA:mature[c]", "RRSB-RRNA:mature[c]", "RRSC-RRNA:mature[c]", "RRSD-RRNA:mature[c]", "RRSE-RRNA:mature[c]", "RRSG-RRNA:mature[c]", "RRSH-RRNA:mature[c]",
	"RRFA-RRNA:mature[c]", "RRFB-RRNA:mature[c]", "RRFC-RRNA:mature[c]", "RRFD-RRNA:mature[c]", "RRFE-RRNA:mature[c]", "RRFF-RRNA:mature[c]", "RRFG-RRNA:mature[c]", "RRFH-RRNA:mature[c]"
	]

_rib23S_IDs = [
	"RRLA-RRNA:mature", "RRLB-RRNA:mature", "RRLC-RRNA:mature","RRLD-RRNA:mature",
	"RRLE-RRNA:mature", "RRLG-RRNA:mature", "RRLH-RRNA:mature"
	]

_rib16S_IDs = [
	"RRSA-RRNA:mature", "RRSB-RRNA:mature", "RRSC-RRNA:mature", "RRSD-RRNA:mature",
	"RRSE-RRNA:mature", "RRSG-RRNA:mature", "RRSH-RRNA:mature"
	]

_rib5S_IDs = [
	"RRFB-RRNA:mature", "RRFC-RRNA:mature", "RRFD-RRNA:mature", "RRFE-RRNA:mature",
	"RRFF-RRNA:mature", "RRFG-RRNA:mature", "RRFH-RRNA:mature"
	]
