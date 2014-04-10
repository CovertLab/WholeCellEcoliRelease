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

		# Constants
		self.elngRate = None
		self.proteinAaCounts = None	# Protein amino acid counts [AA x protein] <-- TODO: Check this
		self.proteinLens = None		# Protein lengths

		super(Translation, self).__init__()

	# Construct object graph, calculate constants
	def initialize(self, sim, kb, kb2):
		super(Translation, self).initialize(sim, kb, kb2)

		# Load parameters
		self.elngRate = kb2.ribosomeElongationRate.to('amino_acid / s').magnitude

		mrnaIDs = kb2.monomerData['rnaId']
		proteinIDs = kb2.monomerData['id']

		# Metabolites
		self.n_aas = len(aaIDs)

		# mRNA, protein monomers
		self.proteinAaCounts = kb2.monomerData['aaCounts'] # TODO: confirm the AA ordering is consistent w/ that used within the process
		self.proteinLens = kb2.monomerData['length']

		# Views
		self.atp = self.bulkMoleculeView('ATP[c]')
		self.adp = self.bulkMoleculeView('ADP[c]')
		self.pi = self.bulkMoleculeView('PI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')

		self.aas = self.bulkMoleculesView(aaIDs)

		# self.mrnas = self.bulkMoleculesView(mrnaIDs)

		self.proteins = self.bulkMoleculesView(proteinIDs)

		self.enzymes = self.bulkMoleculesView(enzIDs)
		self.ribosome23S = self.bulkMoleculesView(rib23S_IDs)
		self.ribosome16S = self.bulkMoleculesView(rib16S_IDs)
		self.ribosome5S = self.bulkMoleculesView(rib5S_IDs)


	def calculateRequest(self):
		nRibosomes = np.min([
			self.ribosome23S.total().sum(),
			self.ribosome16S.total().sum(),
			self.ribosome5S.total().sum(),
			])

		nElongationReactions = np.min([
			nRibosomes * self.elngRate * self.timeStepSec,
			self.aas.total().sum()
			])

		self.aas.requestIs(np.fmin(
			self.aas.total(),
			nElongationReactions//self.n_aas
			))

		self.enzymes.requestAll()
		self.mrnas.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		ribs = np.min([
			self.ribosome23S.total().sum(),
			self.ribosome16S.total().sum(),
			self.ribosome5S.total().sum(),
			])

		# Total synthesis rate
		proteinSynthProb = (self.mrnas.counts() / self.mrnas.counts().sum()).flatten()
		totRate = 1 / np.dot(self.proteinLens, proteinSynthProb) * np.min([
			np.sum(self.aas.counts()),
			ribs * self.elngRate * self.timeStepSec
			])

		# print "Translation totRate: %0.3f" % (totRate)
		newProts = 0
		aasUsed = np.zeros(21)

		proteinsCreated = np.zeros_like(self.proteins.counts())

		# Gillespie-like algorithm
		t = 0
		while True:
			# Choose time step
			t -= np.log(self.randStream.rand()) / totRate
			if t > self.timeStepSec:
				break

			# Check if sufficient metabolic resources to make protein
			newIdx = np.where(self.randStream.mnrnd(1, proteinSynthProb))[0]
			if np.any(self.proteinAaCounts[newIdx, :] > self.aas.counts()):
				break

			# Update metabolites
			self.aas.countsDec(self.proteinAaCounts[newIdx, :].reshape(-1))

			# Increment protein monomer
			proteinsCreated[newIdx] += 1
			newProts += 1
			aasUsed += self.proteinAaCounts[newIdx, :].reshape(-1)

		self.aasUsed = aasUsed

		self.proteins.countsInc(proteinsCreated)


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
