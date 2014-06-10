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
from wholecell.utils.polymerize_new import buildSequences, polymerize, computeMassIncrease, PAD_VALUE


class UniquePolypeptideElongation(wholecell.processes.process.Process):
	""" UniquePolypeptideElongation """

	_name = "UniquePolypeptideElongation"

	# Constructor
	def __init__(self):
		# Parameters
		self.elngRate = None
		self.proteinLengths = None
		self.proteinSequences = None
		self.h2oWeight = None
		self.aaWeightsIncorporated = None

		# Views
		self.activeRibosomes = None
		self.bulkMonomers = None
		self.aas = None
		self.h2o = None
		self.ribosomeSubunits = None

		# Cached values
		self._sequences = None
		self._proteinIndexes = None
		self._peptideLengths = None

		super(UniquePolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniquePolypeptideElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = float(kb.ribosomeElongationRate.to('amino_acid / s').magnitude)

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		proteinIds = kb.monomerData['id']

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

		aaWeights = np.array([
			kb.bulkMolecules[
				kb.bulkMolecules["moleculeId"] == x
				]["mass"].to("fg / mole").magnitude /
			kb.nAvogadro.to("1 / mole").magnitude
			for x in kb.aaIDs
			if len(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == x]["mass"])
			]).flatten()

		self.aaWeightsIncorporated = aaWeights - self.h2oWeight

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		self.aas = self.bulkMoleculesView(aaIds)
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.ribosomeSubunits = self.bulkMoleculesView(enzIds)


	def _buildSequences(self, proteinIndexes, peptideLengths):
		# TODO: cythonize

		# Cache the sequences array used for polymerize, rebuilding if neccesary
		if self._sequences is None or np.any(self._peptideLengths != peptideLengths) or np.any(self._proteinIndexes != proteinIndexes):

			self._sequences = np.empty((proteinIndexes.size, np.int64(self.elngRate)), np.int64)

			for i, (proteinIndex, peptideLength) in enumerate(izip(proteinIndexes, peptideLengths)):
				self._sequences[i, :] = self.proteinSequences[proteinIndex, peptideLength:np.int64(peptideLength + self.elngRate)]

			self._proteinIndexes = proteinIndexes.copy()
			self._peptideLengths = peptideLengths.copy()

			# TODO: if the arrays don't match, try only recomputing the new values
			# and culling the missing entries (this will become important if/when
			# active ribosomes are requested by another process)
		
		return self._sequences


	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		proteinIndexes, peptideLengths = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength'
			)
		
		# sequences = self._buildSequences(proteinIndexes, peptideLengths)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		self.aas.requestIs(
			np.bincount(sequences[sequences != PAD_VALUE])
			)


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

		# sequences = self._buildSequences(proteinIndexes, peptideLengths)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		# Calculate update

		reactionLimit = aaCounts.sum() # TODO: account for energy

		sequenceElongations, aasUsed, nElongations = polymerize(
			sequences,
			aaCounts,
			reactionLimit,
			self.randomState
			)

		# massIncreaseProtein = np.empty_like(massDiffProtein)

		# for i, elongation in enumerate(sequenceElongations):
		# 	massIncreaseProtein[i] = self.aaWeightsIncorporated[sequences[i, :elongation]].sum()

		massIncreaseProtein = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.aaWeightsIncorporated
			)

		updatedMass = massDiffProtein + massIncreaseProtein

		updatedLengths = peptideLengths + sequenceElongations

		didInitialize = (
			(sequenceElongations > 1) &
			(peptideLengths == 0)
			)

		# Update active ribosomes, terminating if neccessary

		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiffProtein = updatedMass
			)

		terminalLengths = self.proteinLengths[proteinIndexes]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedProteins = np.bincount(
			proteinIndexes[didTerminate],
			minlength = self.proteinSequences.shape[0]
			)

		activeRibosomes.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		# Update bulk molecules

		self.aas.countsDec(aasUsed)

		self.bulkMonomers.countsIs(terminatedProteins)

		self.ribosomeSubunits.countsInc(nTerminated)

		self.h2o.countInc(nElongations)

		# Report stalling

		expectedElongations = np.fmin(
			self.elngRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		self.writeToListener("RibosomeStalling", "ribosomeStalls", ribosomeStalls)
