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
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease, PAD_VALUE

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
		self.gtpPerElongation = None

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

		self.aaWeightsIncorporated = kb.translationMonomerWeights

		# TODO: account for ends of peptides in computing weights

		self.gtpPerElongation = kb.gtpPerTranslation

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		self.aas = self.bulkMoleculesView(kb.aaIDs)
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gmp = self.bulkMoleculeView("GMP[c]")
		self.ppi = self.bulkMoleculeView("PPI[c]")

		self.ribosomeSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		proteinIndexes, peptideLengths = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength'
			)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		sequenceHasAA = (sequences != PAD_VALUE)

		self.aas.requestIs(
			np.bincount(sequences[sequenceHasAA])
			)

		self.gtp.requestIs(self.gtpPerElongation * np.fmin(
			sequenceHasAA.sum(),
			self.aas.total().sum()
			))

		# TODO: request water for GTP hydrolysis


	# Calculate temporal evolution
	def evolveState(self):
		aaCounts = self.aas.counts()

		activeRibosomes = self.activeRibosomes.molecules()

		if len(activeRibosomes) == 0:
			return

		proteinIndexes, peptideLengths, massDiffProtein = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength', 'massDiff_protein'
			)

		# Build sequence array

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		# Calculate update

		reactionLimit = self.gtp.count() // self.gtpPerElongation

		sequenceElongations, aasUsed, nElongations = polymerize(
			sequences,
			aaCounts,
			reactionLimit,
			self.randomState
			)

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
			massDiff_protein = updatedMass
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

		gtpUsed = nElongations * self.gtpPerElongation
		self.gtp.countDec(gtpUsed)
		self.gmp.countInc(gtpUsed)
		self.ppi.countInc(gtpUsed)

		# Report stalling

		expectedElongations = np.fmin(
			self.elngRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		self.writeToListener("RibosomeStalling", "ribosomeStalls", ribosomeStalls)
