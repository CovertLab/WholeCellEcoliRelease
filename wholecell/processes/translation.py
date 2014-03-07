#!/usr/bin/env python

"""
Translation

Translation sub-model. Encodes molecular simulation of macromolecular bacterial translation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

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
		self.bulkMoleculesPartition = None

		# Constants
		self.elngRate = 16			# AA/s # TOKB
		self.proteinAaCounts = None	# Protein amino acid counts [AA x protein] <-- TODO: Check this
		self.proteinLens = None		# Protein lengths

		super(Translation, self).__init__()

	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		super(Translation, self).initialize(sim, kb)

		mc = sim.states["BulkMolecules"]

		mrnas = [x for x in kb.rnas if x["monomerId"] != None]
		monomers = [x for x in kb.proteins if len(x["composition"]) == 0 and x["unmodifiedForm"] == None]

		mrnaIDs = [x['id'] + '[c]' for x in mrnas]
		proteinIDs = [x["monomerId"] + "[c]" for x in mrnas]

		# Metabolites
		self.bulkMoleculesPartition.initialize(metIDs + mrnaIDs + proteinIDs + enzIDs)

		self.bulkMoleculesPartition.aas = self.bulkMoleculesPartition.countsView(aaIDs)

		self.n_aas = len(aaIDs)

		self.bulkMoleculesPartition.aasNotSec = self.bulkMoleculesPartition.countsView(aaNotSecIDs)

		self.aasView = mc.countsView(aaIDs)

		self.bulkMoleculesPartition.atpMol = self.bulkMoleculesPartition.countView('ATP[c]')
		self.bulkMoleculesPartition.adpMol = self.bulkMoleculesPartition.countView('ADP[c]')
		self.bulkMoleculesPartition.piMol = self.bulkMoleculesPartition.countView('PI[c]')
		self.bulkMoleculesPartition.h2oMol = self.bulkMoleculesPartition.countView('H2O[c]')
		self.bulkMoleculesPartition.hMol = self.bulkMoleculesPartition.countView('H[c]')

		# mRNA, protein monomers
		self.bulkMoleculesPartition.mrnas = self.bulkMoleculesPartition.countsView(mrnaIDs)
		self.bulkMoleculesPartition.proteins = self.bulkMoleculesPartition.countsView(proteinIDs)
		self.proteinAaCounts = np.array([x["aaCount"] for x in monomers])
		self.proteinLens = np.sum(self.proteinAaCounts, axis = 1)

		# Enzymes
		# TODO: We really want all the associated riboproteins as well (ie we want the complexes)
		self.bulkMoleculesPartition.enzymes = self.bulkMoleculesPartition.countsView(enzIDs)

		self.bulkMoleculesPartition.ribosome23S = self.bulkMoleculesPartition.countsView(rib23S_IDs)
		self.bulkMoleculesPartition.ribosome16S = self.bulkMoleculesPartition.countsView(rib16S_IDs)
		self.bulkMoleculesPartition.ribosome5S = self.bulkMoleculesPartition.countsView(rib5S_IDs)

		self.ribosome23SView = mc.countsView(rib23S_IDs)
		self.ribosome16SView = mc.countsView(rib16S_IDs)
		self.ribosome5SView = mc.countsView(rib5S_IDs)


	def requestBulkMolecules(self):
		self.bulkMoleculesPartition.mrnas.countsIs(1)

		ribs = calcRibosomes(
			self.ribosome23SView.counts(),
			self.ribosome16SView.counts(),
			self.ribosome5SView.counts()
			)

		elng = np.min([
			ribs * self.elngRate * self.timeStepSec,
			np.sum(self.aasView.counts())
			])

		self.bulkMoleculesPartition.aas.countsIs(elng / self.n_aas)


		self.bulkMoleculesPartition.enzymes.countsIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		if not np.any(self.bulkMoleculesPartition.counts()):
			return

		ribs = calcRibosomes(
			self.bulkMoleculesPartition.ribosome23S.counts(),
			self.bulkMoleculesPartition.ribosome16S.counts(),
			self.bulkMoleculesPartition.ribosome5S.counts()
			)

		# print "nRibosomes: %d" % int(calcRibosomes(self.enzyme.counts))
		# Total synthesis rate
		proteinSynthProb = (self.bulkMoleculesPartition.mrnas.counts() / self.bulkMoleculesPartition.mrnas.counts().sum()).flatten()
		totRate = 1 / np.dot(self.proteinLens, proteinSynthProb) * np.min([					# Normalize by average protein length
			np.sum(self.bulkMoleculesPartition.aas.counts()),								# Amino acid limitation
			# np.sum(self.metabolite.counts[self.metabolite.idx["atp"]]) / 2,							# GTP (energy) limitation
			ribs * self.elngRate * self.timeStepSec
			# self.enzyme.counts[self.enzyme.idx["ribosome70S"]] * self.elngRate * self.timeStepSec		# Ribosome capacity
			])

		# print "Translation totRate: %0.3f" % (totRate)
		newProts = 0
		aasUsed = np.zeros(21)

		proteinsCreated = np.zeros_like(self.bulkMoleculesPartition.proteins.counts())

		# Gillespie-like algorithm
		t = 0
		while True:
			# Choose time step
			t -= np.log(self.randStream.rand()) / totRate
			if t > self.timeStepSec:
				break

			# Check if sufficient metabolic resources to make protein
			newIdx = np.where(self.randStream.mnrnd(1, proteinSynthProb))[0]
			if \
				np.any(self.proteinAaCounts[newIdx, :] > self.bulkMoleculesPartition.aas.counts()):
				# np.any(self.proteinAaCounts[newIdx, :] > self.metabolite.counts[self.metabolite.idx["aas"]]) or \
				# 2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["atp"]] or \
				# 2 * self.proteinLens[newIdx] > self.metabolite.counts[self.metabolite.idx["h2o"]]:
					break

			# Update metabolites
			#self.metabolite.counts[self.metabolite.idx["aas"]] -= self.proteinAaCounts[newIdx, :].reshape(-1)
			self.bulkMoleculesPartition.aas.countsDec(self.proteinAaCounts[newIdx, :].reshape(-1))
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
		# print "Translation numActiveRibs (total): %d (%d)" % (int(np.sum(aasUsed) / self.elngRate / self.timeStepSec), int(calcRibosomes(self.enzyme.counts)))

		self.bulkMoleculesPartition.proteins.countsInc(proteinsCreated)


def calcRibosomes(counts23S, counts16S, counts5S):
	return np.min([
		np.sum(counts23S),
		np.sum(counts16S),
		np.sum(counts5S)
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
