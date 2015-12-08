#!/usr/bin/env python

"""
PolypeptideInitiation

Polypeptide initiation sub-model.

TODO:
- translate off of transcribed transcription units
- co-transcriptional translation
- retain mRNA while translating
- match observed active ribosome fraction
- use actual ribosomes (all isoforms, all subunits) instead of representative rRNAs

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.fitting import normalize

import itertools

class PolypeptideInitiation(wholecell.processes.process.Process):
	""" PolypeptideInitiation """

	_name = "PolypeptideInitiation"

	# Constructor
	def __init__(self):
		# Parameters
		self.proteinLens = None

		# Views

		self.activeRibosomes = None
		self.ribosome30S = None
		self.ribosome50S = None
		self.mRnas = None

		super(PolypeptideInitiation, self).__init__()


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(PolypeptideInitiation, self).initialize(sim, sim_data)

		# Load parameters

		mrnaIds = sim_data.process.translation.monomerData["rnaId"]

		self.proteinLens = sim_data.process.translation.monomerData["length"].asNumber()

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeGroups.s30_fullComplex[0])
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeGroups.s50_fullComplex[0])

		self.mRnas = self.bulkMoleculesView(mrnaIds)
		self.translationEfficiencies = normalize(sim_data.process.translation.translationEfficienciesByMonomer)


	def calculateRequest(self):
		self.ribosome30S.requestAll()
		self.ribosome50S.requestAll()

		self.mRnas.requestAll()

	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRibosomeCount = np.min([
			self.ribosome30S.count().sum(),
			self.ribosome50S.count().sum(),
			])

		if inactiveRibosomeCount == 0:
			return

		proteinInitProb = normalize(
			self.mRnas.counts() * self.translationEfficiencies
			)

		nNewProteins = self.randomState.multinomial(
			inactiveRibosomeCount,
			proteinInitProb
			)

		nonzeroCount = (nNewProteins > 0)

		assert nNewProteins.sum() == inactiveRibosomeCount

		# Build list of protein indexes

		proteinIndexes = np.empty(inactiveRibosomeCount, np.int64)

		startIndex = 0
		for proteinIndex, counts in itertools.izip(
				np.arange(nNewProteins.size)[nonzeroCount],
				nNewProteins[nonzeroCount],
				):

			proteinIndexes[startIndex:startIndex+counts] = proteinIndex

			startIndex += counts

		# Create the active ribosomes

		activeRibosomes = self.activeRibosomes.moleculesNew(
			"activeRibosome",
			inactiveRibosomeCount
			)

		activeRibosomes.attrIs(
			proteinIndex = proteinIndexes
			)

		self.ribosome30S.countDec(nNewProteins.sum())
		self.ribosome50S.countDec(nNewProteins.sum())

		self.writeToListener("RibosomeData", "didInitialize", nNewProteins.sum())
