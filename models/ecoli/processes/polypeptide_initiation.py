#!/usr/bin/env python

"""
PolypeptideInitiation

Polypeptide initiation sub-model.

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

	def __init__(self):
		super(PolypeptideInitiation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideInitiation, self).initialize(sim, sim_data)

		# Load parameters
		mrnaIds = sim_data.process.translation.monomerData["rnaId"]
		self.proteinLens = sim_data.process.translation.monomerData["length"].asNumber()
		self.translationEfficiencies = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

		# Create view on to active 70S ribosomes
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		# Create views onto bulk 30S and 50S ribosomal subunits
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeGroups.s30_fullComplex[0])
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeGroups.s50_fullComplex[0])

		# Create view onto bulk mRNAs
		self.mRnas = self.bulkMoleculesView(mrnaIds)

	def calculateRequest(self):
		self.ribosome30S.requestAll()
		self.ribosome50S.requestAll()
		self.mRnas.requestAll()

	def evolveState(self):
		# Calculate number of ribosomes that could potentially be initalized based on
		# counts of free 30S and 50S subunits
		inactiveRibosomeCount = np.min([
			self.ribosome30S.count().sum(),
			self.ribosome50S.count().sum(),
			])

		if inactiveRibosomeCount == 0:
			return

		# Calculate initiation probabilities for ribosomes based on mRNA counts and associated
		# mRNA translational efficiencies
		proteinInitProb = normalize(
			self.mRnas.counts() * self.translationEfficiencies
			)

		# Sample multinomial distribution to determine which mRNAs have full 70S
		# ribosomes initalized on them
		nNewProteins = self.randomState.multinomial(
			inactiveRibosomeCount,
			proteinInitProb
			)

		# Check that sampling produced expected result
		assert nNewProteins.sum() == inactiveRibosomeCount

		# Each ribosome is assigned a protein index for the protein that corresponds to the
		# polypeptide it will polymerize. This is done in blocks of protein ids for efficiency.
		proteinIndexes = np.empty(inactiveRibosomeCount, np.int64)
		nonzeroCount = (nNewProteins > 0)
		startIndex = 0
		for proteinIndex, counts in itertools.izip(
				np.arange(nNewProteins.size)[nonzeroCount],
				nNewProteins[nonzeroCount],
				):

			proteinIndexes[startIndex:startIndex+counts] = proteinIndex

			startIndex += counts

		# Create active 70S ribosomes and assign their protein indexes calculated above
		activeRibosomes = self.activeRibosomes.moleculesNew(
			"activeRibosome",
			inactiveRibosomeCount
			)

		activeRibosomes.attrIs(
			proteinIndex = proteinIndexes
			)

		# Decrement free 30S and 70S ribosomal subunit counts
		self.ribosome30S.countDec(nNewProteins.sum())
		self.ribosome50S.countDec(nNewProteins.sum())

		# Write number of initalized ribosomes to listener
		self.writeToListener("RibosomeData", "didInitialize", nNewProteins.sum())
