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
		
		# Partitions
		self.metabolitePartition = None
		self.mrnaPartition = None
		self.proteinPartition = None
		self.enzymePartition = None

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
		self.metabolitePartition = mc.addPartition(self, _metIDs, self.calcReqMetabolites)

		self.metabolitePartition.aas = self.metabolitePartition.countsBulkViewNew(_aaIDs)

		self.n_aas = len(_aaIDs)

		self.metabolitePartition.aasNotSec = self.metabolitePartition.countsBulkViewNew(_aaNotSecIDs)

		self.aas = mc.countsBulkViewNew(_aaIDs)

		self.metabolitePartition.atpMol = self.metabolitePartition.molecule('ATP', 'merged')
		self.metabolitePartition.adpMol = self.metabolitePartition.molecule('ADP', 'merged')
		self.metabolitePartition.piMol = self.metabolitePartition.molecule('PI', 'merged')
		self.metabolitePartition.h2oMol = self.metabolitePartition.molecule('H2O', 'merged')
		self.metabolitePartition.hMol = self.metabolitePartition.molecule('H', 'merged')

		# mRNA, protein monomer
		mrnas = [x for x in kb.rnas if x["monomerId"] != None]
		monomers = [x for x in kb.proteins if len(x["composition"]) == 0 and x["unmodifiedForm"] == None]
		self.mrnaPartition = mc.addPartition(self,[x["id"] + "[c]" for x in mrnas], self.calcReqMrna)
		self.proteinPartition = mc.addPartition(self,[x["monomerId"] + ":nascent[c]" for x in mrnas], self.calcReqProtein)
		self.proteinAaCounts = numpy.array([x["aaCount"] for x in monomers])
		self.proteinLens = numpy.sum(self.proteinAaCounts, axis = 1)

		# Enzymes
		# TODO: We really want all the associated riboproteins as well (ie we want the complexes)
		self.enzymePartition = mc.addPartition(self, _enzIDs, self.calcReqEnzyme)

		self.enzymePartition.ribosome23S = self.enzymePartition.countsBulkViewNew(_rib23S_IDs)

		self.enzymePartition.ribosome16S = self.enzymePartition.countsBulkViewNew(_rib16S_IDs)

		self.enzymePartition.ribosome5S = self.enzymePartition.countsBulkViewNew(_rib5S_IDs)

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
		if not numpy.any(self.mrnaPartition.countsBulk()):
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
	"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
	"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	]

_aaNotSecIDs = [
	"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",
	"LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	]

_enzIDs = [
	"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]", "RRLD-RRNA[c]", "RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]",
	"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]", "RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]",
	"RRFA-RRNA[c]", "RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]", "RRFE-RRNA[c]", "RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
	]

_rib23S_IDs = [
	"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]","RRLD-RRNA[c]",
	"RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]"
	]

_rib16S_IDs = [
	"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]",
	"RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]"
	]

_rib5S_IDs = [
	"RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]", "RRFE-RRNA[c]",
	"RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
	]
