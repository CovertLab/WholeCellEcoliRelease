#!/usr/bin/env python

"""
UniquePolypeptideElongation

Translation elongation sub-model.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize_new import polymerize, PAD_VALUE


class UniquePolypeptideElongation(wholecell.processes.process.Process):
	""" UniquePolypeptideElongation """

	_name = "UniquePolypeptideElongation"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None
		self.proteinIds = None

		super(UniquePolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniquePolypeptideElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = float(kb.ribosomeElongationRate.to('amino_acid / s').magnitude)

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		self.proteinIds = kb.monomerData['id']

		self.proteinLengths = kb.monomerData["length"].magnitude

		self.proteinSequences = kb.translationSequences

		aaIds = kb.aaIDs[:]

		# TODO: Remove hack of deleting selenocysteine this way
		selenocysteineIdx = aaIds.index("SEC-L[c]")
		del aaIds[selenocysteineIdx]

		# # TODO: refactor mass updates

		self.h2oWeight = (
			kb.bulkMolecules[
				kb.bulkMolecules["moleculeId"] == "H2O[c]"
				]["mass"].to("fg / mole").magnitude /
			kb.nAvogadro.to("1 / mole").magnitude
			)

		self.aaWeights = np.array([
			kb.bulkMolecules[
				kb.bulkMolecules["moleculeId"] == x
				]["mass"].to("fg / mole").magnitude /
			kb.nAvogadro.to("1 / mole").magnitude
			for x in kb.aaIDs
			if len(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == x]["mass"])
			])

		self.aaWeightsIncorporated = self.aaWeights - self.h2oWeight

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(self.proteinIds)

		self.aas = self.bulkMoleculesView(aaIds)
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.ribosomeSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		proteinIndexes, peptideLengths = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength'
			)
		
		# HACK DO NOT COMMIT
		self.sequences = np.empty((proteinIndexes.size, np.int64(self.elngRate)), np.int64)

		for i, (proteinIndex, peptideLength) in enumerate(izip(proteinIndexes, peptideLengths)):
			self.sequences[i, :] = self.proteinSequences[proteinIndex, peptideLength:np.int64(peptideLength + self.elngRate)]

		self.aas.requestIs(
			np.bincount(self.sequences[self.sequences != PAD_VALUE])
			)

		# self.aas.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		aaCounts = self.aas.counts()

		activeRibosomes = self.activeRibosomes.molecules()

		if len(activeRibosomes) == 0:
			return

		proteinIndexes, peptideLengths, massDiffProtein = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength', 'massDiffProtein'
			)

		# Build sequence array

		# HACK DO NOT COMMIT
		# sequences = np.empty((proteinIndexes.size, np.int64(self.elngRate)), np.int64)

		# for i, (proteinIndex, peptideLength) in enumerate(izip(proteinIndexes, peptideLengths)):
		# 	sequences[i, :] = self.proteinSequences[proteinIndex, peptideLength:np.int64(peptideLength + self.elngRate)]

		# Calculate update

		reactionLimit = aaCounts.sum() # TODO: account for energy

		sequenceElongation, aasUsed, nElongations = polymerize(
			self.sequences,
			aaCounts,
			reactionLimit,
			self.randStream
			)

		updatedMass = massDiffProtein + np.array([
			self.aaWeightsIncorporated[self.sequences[i, :elongation]].sum()
			for i, elongation in enumerate(sequenceElongation)
			])

		updatedLengths = peptideLengths + sequenceElongation

		didInitialize = (
			(sequenceElongation > 1) &
			(peptideLengths == 0)
			)

		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiffProtein = updatedMass
			)

		terminatedProteins = np.zeros_like(self.bulkMonomers.counts())

		didTerminate = (updatedLengths == self.proteinLengths[proteinIndexes])

		for moleculeIndex, molecule in enumerate(activeRibosomes):
			if didTerminate[moleculeIndex]:
				terminatedProteins[molecule.attr('proteinIndex')] += 1
				self.activeRibosomes.moleculeDel(molecule)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		self.aas.countsDec(aasUsed)

		self.bulkMonomers.countsIs(terminatedProteins)

		self.ribosomeSubunits.countsInc(nTerminated)

		self.h2o.countInc(nElongations)


def getWorkAssignment(dataSize, thisTask, totalTasks):
	startPos = sum(dataSize // totalTasks + (i < (dataSize % totalTasks)) for i in xrange(thisTask))
	nElements = dataSize // totalTasks + (thisTask < (dataSize % totalTasks))

	return startPos, nElements