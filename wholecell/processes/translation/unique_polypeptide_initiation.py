#!/usr/bin/env python

"""
UniquePolypeptideInitiation

Polypeptide initiation sub-model.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

import itertools

class UniquePolypeptideInitiation(wholecell.processes.process.Process):
	""" UniquePolypeptideInitiation """

	_name = "UniquePolypeptideInitiation"

	# Constructor
	def __init__(self):
		# Constants

		self.proteinAACounts = None

		super(UniquePolypeptideInitiation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniquePolypeptideInitiation, self).initialize(sim, kb)

		# Load parameters

		mrnaIds = kb.monomerData["rnaId"]
		proteinIds = kb.monomerData["id"]
		aaIds = kb.aaIDs[:]

		# TODO: Remove hack of deleting selenocysteine this way
		selenocysteineIdx = aaIds.index("SEC-L[c]")
		del aaIds[selenocysteineIdx]

		self.proteinAACounts = np.delete(
			kb.monomerData["aaCounts"].magnitude, selenocysteineIdx, 1
			)
		self.proteinLens = kb.monomerData["length"].magnitude

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		self.rRna23S = self.bulkMoleculesView([rib23S_IDs[0]])
		self.rRna16S = self.bulkMoleculesView([rib16S_IDs[0]])
		self.rRna5S = self.bulkMoleculesView([rib5S_IDs[0]])

		self.mRnas = self.bulkMoleculesView(mrnaIds)


	def calculateRequest(self):
		self.rRna23S.requestAll()
		self.rRna16S.requestAll()
		self.rRna5S.requestAll()

		self.mRnas.requestAll()

	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRibosomes = np.min([
			self.rRna23S.counts().sum(),
			self.rRna16S.counts().sum(),
			self.rRna5S.counts().sum()
			])

		proteinInitProb = (
			self.mRnas.counts() /
			self.mRnas.counts().sum()
			).flatten()	# TODO: Is this .flatten() necessary?

		nNewProteins = self.randStream.mnrnd(
			inactiveRibosomes,
			proteinInitProb
			)

		nonzeroCount = (nNewProteins > 0)

		for protIdx, nNew, aaCounts in itertools.izip(
			np.arange(nNewProteins.size)[nonzeroCount],
			nNewProteins[nonzeroCount],
			self.proteinAACounts[nonzeroCount]
			):

			self.activeRibosomes.moleculesNew(
				"activeRibosome",
				nNew,
				proteinIndex = protIdx,
				requiredAAs = aaCounts
				)

		self.rRna23S.countsDec(nNewProteins.sum())
		self.rRna16S.countsDec(nNewProteins.sum())
		self.rRna5S.countsDec(nNewProteins.sum())
		print nNewProteins.sum(), inactiveRibosomes

# TODO: To kb
rib23S_IDs = [
	"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]","RRLD-RRNA[c]",
	"RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]"
	]

rib16S_IDs = [
	"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]",
	"RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]"
	]

rib5S_IDs = [
	"RRFA-RRNA[c]", "RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]",
	"RRFE-RRNA[c]", "RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
	]
