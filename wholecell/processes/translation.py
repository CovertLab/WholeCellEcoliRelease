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

	_name = "Translation"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None
		self.proteinAaCounts = None	# Protein amino acid counts [AA x protein] <-- TODO: Check this
		self.proteinLens = None		# Protein lengths

		super(Translation, self).__init__()

	# Construct object graph, calculate constants
	def initialize(self, sim, kb):
		super(Translation, self).initialize(sim, kb)

		# Load parameters
		self.elngRate = kb.ribosomeElongationRate.to('amino_acid / s').magnitude

		mrnaIDs = kb.monomerData['rnaId']
		proteinIDs = kb.monomerData['id']

		# Metabolites
		self.n_aas = len(aaIDs)

		oneToThreeMapping = dict((
			("A", "ALA-L[c]"), ("R", "ARG-L[c]"), ("N", "ASN-L[c]"), ("D", "ASP-L[c]"),
			("C", "CYS-L[c]"), ("E", "GLU-L[c]"), ("Q", "GLN-L[c]"), ("G", "GLY[c]"),
			("H", "HIS-L[c]"), ("I", "ILE-L[c]"), ("L", "LEU-L[c]"), ("K", "LYS-L[c]"),
			("M", "MET-L[c]"), ("F", "PHE-L[c]"), ("P", "PRO-L[c]"), ("S", "SER-L[c]"),
			("T", "THR-L[c]"), ("W", "TRP-L[c]"), ("Y", "TYR-L[c]"), ("U", "SEC-L[c]"),
			("V", "VAL-L[c]")
		)) # TOKB

		# mRNA, protein monomers
		self.proteinAaCounts = kb.monomerData['aaCounts'] # TODO: confirm the AA ordering is consistent w/ that used within the process
		self.proteinLens = kb.monomerData['length']
		self.avgProteinLength = np.mean(self.proteinLens)

		self.fracActiveRibosomes = kb.fracActiveRibosomes.magnitude

		# Views
		self.atp = self.bulkMoleculeView('ATP[c]')
		self.adp = self.bulkMoleculeView('ADP[c]')
		self.pi = self.bulkMoleculeView('PI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')

		self.aas = self.bulkMoleculesView(
			[oneToThreeMapping[x] if x != "U" else "H2O[c]" for x in kb._aaWeights.iterkeys()] #
			)

		self.mrnas = self.bulkMoleculesView(mrnaIDs)

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

		# self.aas.requestIs(np.fmin(
		# 	self.aas.total(),
		# 	nElongationReactions//self.n_aas
		# 	))

		self.aas.requestAll()
		self.enzymes.requestAll()
		self.mrnas.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		ribs = np.min([
			self.ribosome23S.total().sum(),
			self.ribosome16S.total().sum(),
			self.ribosome5S.total().sum(),
			])

		proteinSynthProb = (self.mrnas.counts() / self.mrnas.counts().sum()).flatten()

		enzLimit = ribs * self.elngRate * self.timeStepSec

		# nProteinsToCreate = int(enzLimit / self.avgProteinLength)

		# aasShape = self.aas.counts().shape

		# proteinsCreated = np.zeros_like(self.proteins.counts())

		# while enzLimit > 0 and nProteinsToCreate > 0:
		# 	# print nProteinsToCreate
		# 	if not np.any(
		# 		np.all(
		# 			self.aas.counts() > self.proteinAaCounts,
		# 			axis = 1
		# 			)
		# 		):
		# 		break

		# 	if not np.any(enzLimit > np.sum(self.proteinAaCounts, axis = 1)):
		# 		break

		# 	# If the probabilities of being able to synthesize are sufficiently low, exit the loop
		# 	if np.sum(proteinSynthProb[np.all(self.aas.counts() > self.proteinAaCounts, axis = 1)]) < 1e-3:
		# 		break

		# 	if np.sum(proteinSynthProb[enzLimit > np.sum(self.proteinAaCounts, axis = 1)]) < 1e-3:
		# 		break

		# 	newIdxs = np.where(self.randStream.mnrnd(nProteinsToCreate, proteinSynthProb))[0]

		# 	if np.any(self.aas.counts() < np.sum(self.proteinAaCounts[newIdxs, :], axis = 0)):
		# 		if nProteinsToCreate > 0:
		# 			nProteinsToCreate //= 2
		# 			continue
		# 		else:
		# 			break

		# 	if enzLimit < np.sum(np.sum(self.proteinAaCounts[newIdxs, :], axis = 0)):
		# 		if nProteinsToCreate > 0:
		# 			nProteinsToCreate //= 2
		# 			continue
		# 		else:
		# 			break

		# 	enzLimit -= np.sum(self.proteinAaCounts[newIdxs, :])

		# 	self.aas.countsDec(
		# 		np.sum(self.proteinAaCounts[newIdxs, :], axis = 0).reshape(aasShape)
		# 		)

		# 	# TODO: Handle water, etc
		# 	proteinsCreated[newIdxs] += 1

		from wholecell.utils.simple_polymerize import simplePolymerize

		aaCounts, proteinsCreated = simplePolymerize(
			self.proteinAaCounts,
			enzLimit,
			self.aas.counts(),
			proteinSynthProb,
			self.randStream
			)

		self.aas.countsIs(aaCounts)

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
