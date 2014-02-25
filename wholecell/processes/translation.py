#!/usr/bin/env python

"""
Translation

Translation sub-model. Encodes molecular simulation of macromolecular bacterial translation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.processes.process

class Translation(wholecell.processes.process.Process):
	""" Translation """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "Translation",
			"name": "Translation"
		}
		
		# Partitions
		self.mcPartition = None

		# Constants
		self.elngRate = 16			# AA/s # TOKB
		self.proteinAaCounts = None	# Protein amino acid counts [AA x protein] <-- TODO: Check this
		self.proteinLens = None		# Protein lengths

		super(Translation, self).__init__()

	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		super(Translation, self).initialize(sim, kb)

		mc = sim.states["MoleculeCounts"]

		mrnas = [x for x in kb.rnas if x["monomerId"] != None]
		monomers = [x for x in kb.proteins if len(x["composition"]) == 0 and x["unmodifiedForm"] == None]

		mrnaIDs = [x['id'] + '[c]' for x in mrnas]
		proteinIDs = [x["monomerId"] + ":nascent[c]" for x in mrnas]

		# Metabolites
		self.mcPartition.initialize(metIDs + mrnaIDs + proteinIDs + enzIDs)

		self.mcPartition.aas = self.mcPartition.countsBulkViewNew(aaIDs)

		self.n_aas = len(aaIDs)

		self.mcPartition.aasNotSec = self.mcPartition.countsBulkViewNew(aaNotSecIDs)

		self.aasView = mc.countsBulkViewNew(aaIDs)

		self.mcPartition.atpMol = self.mcPartition.molecule('ATP[c]')
		self.mcPartition.adpMol = self.mcPartition.molecule('ADP[c]')
		self.mcPartition.piMol = self.mcPartition.molecule('PI[c]')
		self.mcPartition.h2oMol = self.mcPartition.molecule('H2O[c]')
		self.mcPartition.hMol = self.mcPartition.molecule('H[c]')

		# mRNA, protein monomers
		self.mcPartition.mrnas = self.mcPartition.countsBulkViewNew(mrnaIDs)
		self.mcPartition.proteins = self.mcPartition.countsBulkViewNew(proteinIDs)
		self.proteinAaCounts = numpy.array([x["aaCount"] for x in monomers])
		self.proteinLens = numpy.sum(self.proteinAaCounts, axis = 1)

		# Enzymes
		# TODO: We really want all the associated riboproteins as well (ie we want the complexes)
		self.mcPartition.enzymes = self.mcPartition.countsBulkViewNew(enzIDs)

		self.mcPartition.ribosome23S = self.mcPartition.countsBulkViewNew(rib23S_IDs)
		self.mcPartition.ribosome16S = self.mcPartition.countsBulkViewNew(rib16S_IDs)
		self.mcPartition.ribosome5S = self.mcPartition.countsBulkViewNew(rib5S_IDs)

		self.ribosome23SView = mc.countsBulkViewNew(rib23S_IDs)
		self.ribosome16SView = mc.countsBulkViewNew(rib16S_IDs)
		self.ribosome5SView = mc.countsBulkViewNew(rib5S_IDs)


	def requestMoleculeCounts(self):
		self.mcPartition.mrnas.countsBulkIs(1)

		ribs = calcRibosomes(
			self.ribosome23SView.countsBulk(),
			self.ribosome16SView.countsBulk(),
			self.ribosome5SView.countsBulk()
			)

		elng = numpy.min([
			ribs * self.elngRate * self.timeStepSec,
			numpy.sum(self.aasView.countsBulk())
			])

		self.mcPartition.aas.countsBulkIs(elng / self.n_aas)


		self.mcPartition.enzymes.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		if not numpy.any(self.mcPartition.countsBulk()):
			return

		ribs = calcRibosomes(
			self.mcPartition.ribosome23S.countsBulk(),
			self.mcPartition.ribosome16S.countsBulk(),
			self.mcPartition.ribosome5S.countsBulk()
			)

		# print "nRibosomes: %d" % int(calcRibosomes(self.enzyme.counts))
		# Total synthesis rate
		proteinSynthProb = (self.mcPartition.mrnas.countsBulk() / numpy.sum(self.mcPartition.mrnas.countsBulk())).flatten()
		totRate = 1 / numpy.dot(self.proteinLens, proteinSynthProb) * numpy.min([					# Normalize by average protein length
			numpy.sum(self.mcPartition.aas.countsBulk()),								# Amino acid limitation
			# numpy.sum(self.metabolite.counts[self.metabolite.idx["atp"]]) / 2,							# GTP (energy) limitation
			ribs * self.elngRate * self.timeStepSec
			# self.enzyme.counts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec		# Ribosome capacity
			])

		# print "Translation totRate: %0.3f" % (totRate)
		newProts = 0
		aasUsed = numpy.zeros(21)

		proteinsCreated = numpy.zeros_like(self.mcPartition.proteins.countsBulk())

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
				numpy.any(self.proteinAaCounts[newIdx, :] > self.mcPartition.aas.countsBulk()):
				# numpy.any(self.proteinAaCounts[newIdx, :] > self.metabolite.counts[self.metabolite.idx["aas"]]) or \
				# 2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["atp"]] or \
				# 2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["h2o"]]:
					break

			# Update metabolites
			#self.metabolite.counts[self.metabolite.idx["aas"]] -= self.proteinAaCounts[newIdx, :].reshape(-1)
			self.mcPartition.aas.countsBulkDec(self.proteinAaCounts[newIdx, :].reshape(-1))
			# self.metabolite.counts[self.metabolite.idx["h2o"]] += self.proteinLens[newIdx] - 1

			# self.metabolite.counts[self.metabolite.idx["atp"]] -= 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["h2o"]] -= 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["adp"]] += 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["pi"]] += 2 * self.proteinLens[newIdx]
			# self.metabolite.counts[self.metabolite.idx["h"]] += 2 * self.proteinLens[newIdx]

			# Increment protein monomer
			#self.protein.counts[newIdx] += 1
			proteinsCreated[newIdx] += 1
			newProts += 1
			aasUsed += self.proteinAaCounts[newIdx, :].reshape(-1)

		self.aasUsed = aasUsed
		# print "Translation newProts: %d" % newProts
		# print "Translation aasUsed: %s" % str(aasUsed)
		# print "Translation numActiveRibs (total): %d (%d)" % (int(numpy.sum(aasUsed) / self.elngRate / self.timeStepSec), int(calcRibosomes(self.enzyme.counts)))

		self.mcPartition.proteins.countsBulkInc(proteinsCreated)


def calcRibosomes(counts23S, counts16S, counts5S):
	return numpy.min([
		numpy.sum(counts23S),
		numpy.sum(counts16S),
		numpy.sum(counts5S)
		])


metIDs = ["ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]",
	"GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
	"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]",
	"THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	#"FMET[c]", # TODO: Re-add
	#"GTP[c]", "GDP[c]", "PI[c]",  "H2O[c]", "H[c]"
	"ATP[c]", "ADP[c]", "PI[c]",  "H2O[c]", "H[c]"]

aaIDs = [
	"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
	"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	]

aaNotSecIDs = [
	"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",
	"LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	]

enzIDs = [
	"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]", "RRLD-RRNA[c]", "RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]",
	"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]", "RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]",
	"RRFA-RRNA[c]", "RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]", "RRFE-RRNA[c]", "RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
	]

rib23S_IDs = [
	"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]","RRLD-RRNA[c]",
	"RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]"
	]

rib16S_IDs = [
	"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]",
	"RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]"
	]

rib5S_IDs = [
	"RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]", "RRFE-RRNA[c]",
	"RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
	]
