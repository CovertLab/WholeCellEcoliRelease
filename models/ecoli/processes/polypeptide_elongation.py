#!/usr/bin/env python

"""
PolypeptideElongation

Translation elongation sub-model.

TODO:
- see the initiation process for more TODOs

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease, PAD_VALUE
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		# Parameters
		self.elngRate = None
		self.proteinLengths = None
		self.proteinSequences = None
		self.h2oWeight = None
		self.aaWeightsIncorporated = None
		self.gtpPerElongation = None
		self.synthetase_turnover = None

		# Views
		self.activeRibosomes = None
		self.bulkMonomers = None
		self.aas = None
		self.h2o = None
		self.trna_groups = None
		self.synthetase_groups = None
		self.ribosomeSubunits = None

		super(PolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = float(kb.ribosomeElongationRate.asNumber(units.aa / units.s)) * self.timeStepSec

		self.aa_trna_groups = kb.aa_trna_groups
		self.aa_synthetase_groups = kb.aa_synthetase_groups
		self.synthetase_turnover = kb.trna_synthetase_rates.asNumber(units.aa/units.s)

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		proteinIds = kb.monomerData['id']

		self.proteinLengths = kb.monomerData["length"].asNumber()

		self.proteinSequences = kb.translationSequences

		self.aaWeightsIncorporated = kb.translationMonomerWeights

		self.endWeight = kb.translationEndWeight

		self.gtpPerElongation = kb.gtpPerTranslation

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		self.aas = self.bulkMoleculesView(kb.aaIDs)
		self.trna_groups = [self.bulkMoleculesView(x) for x in self.aa_trna_groups.itervalues()]
		self.synthetase_groups = [self.bulkMoleculesView(x) for x in self.aa_synthetase_groups.itervalues()]
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h   = self.bulkMoleculeView("H[c]")

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

		aasRequested = np.bincount(sequences[sequenceHasAA])

		self.aas.requestIs(
			aasRequested
			)

		# Should essentially request all tRNAs
		# and all synthetases
		trnasRequested = aasRequested
		for i,group in enumerate(self.trna_groups):
			group.requestIs(trnasRequested[i])
		synthetaseRequested = aasRequested
		for i,group in enumerate(self.synthetase_groups):
			group.requestIs(synthetaseRequested[i])

		gtpsHydrolyzed = np.int64(np.ceil(
			self.gtpPerElongation * np.fmin(
				sequenceHasAA.sum(),
				self.aas.total().sum()
				)
			))

		self.gtp.requestIs(gtpsHydrolyzed)

		self.h2o.requestIs(gtpsHydrolyzed) # note: this is roughly a 2x overestimate


	# Calculate temporal evolution
	def evolveState(self):
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

		# Calculate elongation resource capacity

		aaCountInSequence = np.bincount(sequences[(sequences != PAD_VALUE)])
		aaCounts = self.aas.counts()
		trnasCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.trna_groups],dtype = np.int64)
		synthetaseCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.synthetase_groups],dtype = np.int64)
		elongationResourceCapacity = np.minimum(aaCounts, synthetaseCapacity, trnasCapacity)

		# Calculate update

		reactionLimit = self.gtp.count() // self.gtpPerElongation

		sequenceElongations, aasUsed, nElongations = polymerize(
			sequences,
			aaCounts, # elongationResourceCapacity,
			reactionLimit,
			self.randomState
			)

		massIncreaseProtein = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.aaWeightsIncorporated
			)

		updatedLengths = peptideLengths + sequenceElongations

		didInitialize = (
			(sequenceElongations > 1) &
			(peptideLengths == 0)
			)

		updatedMass = massDiffProtein + massIncreaseProtein

		updatedMass[didInitialize] += self.endWeight

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

		self.bulkMonomers.countsInc(terminatedProteins)

		self.ribosomeSubunits.countsInc(nTerminated)

		self.h2o.countInc(nElongations - nInitialized)

		gtpUsed = np.int64(stochasticRound(
			self.randomState, 
			nElongations * self.gtpPerElongation
			))

		self.gtp.countDec(gtpUsed)
		self.gdp.countInc(gtpUsed)
		self.pi.countInc(gtpUsed)
		self.h.countInc(gtpUsed)

		self.h2o.countDec(gtpUsed)

		# Report stalling information

		expectedElongations = np.fmin(
			self.elngRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		self.writeToListener("RibosomeStalling", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeStalling", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeStalling", "aaCounts", aaCounts)
		self.writeToListener("RibosomeStalling", "trnasCapacity", trnasCapacity)
		self.writeToListener("RibosomeStalling", "synthetaseCapacity", synthetaseCapacity)
