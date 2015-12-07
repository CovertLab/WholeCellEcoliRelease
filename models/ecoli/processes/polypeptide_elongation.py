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

SYNTHETASE_KM_SCALE = 0.1

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
		self.ribosome30S = None
		self.ribosome50S = None

		super(PolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = float(kb.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s)) * self.timeStepSec
		self.elngRate = int(round(self.elngRate)) # TODO: Make this less of a hack by implementing in the KB

		self.nAvogadro = kb.constants.nAvogadro
		self.cellDensity = kb.constants.cellDensity

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		proteinIds = kb.process.translation.monomerData['id']

		self.proteinLengths = kb.process.translation.monomerData["length"].asNumber()

		self.proteinSequences = kb.process.translation.translationSequences

		self.aaWeightsIncorporated = kb.process.translation.translationMonomerWeights

		self.endWeight = kb.process.translation.translationEndWeight

		self.gtpPerElongation = kb.constants.gtpPerTranslation

		##########
		aaIdxs = [kb.process.metabolism.metabolitePoolIDs.index(aaID) for aaID in kb.moleculeGroups.aaIDs]
		aaConcentrations = kb.process.metabolism.metabolitePoolConcentrations[aaIdxs]
		total_aa_concentration = units.sum(aaConcentrations)
		self.saturation_km = SYNTHETASE_KM_SCALE * total_aa_concentration
		##########

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		self.aas = self.bulkMoleculesView(kb.moleculeGroups.aaIDs)
		self.h2o = self.bulkMoleculeView('WATER[c]')

		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("Pi[c]")
		self.h   = self.bulkMoleculeView("PROTON[c]")

		self.ribosome30S = self.bulkMoleculeView(kb.moleculeGroups.s30_fullComplex[0])
		self.ribosome50S = self.bulkMoleculeView(kb.moleculeGroups.s50_fullComplex[0])

		###### VARIANT CODE #######
		self.translationSaturation = kb.translationSaturation
		###### VARIANT CODE #######

	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		if len(activeRibosomes) == 0:
			return

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

		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		if self.translationSaturation:
			cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
			cellVolume = cellMass / self.cellDensity
			aaTotalConc = units.sum((1 / self.nAvogadro) * (1 / cellVolume) * self.aas.total())

			translation_machinery_saturation = (aaTotalConc / (self.saturation_km + aaTotalConc))
			translation_machinery_saturation = units.convertNoUnitToNumber(translation_machinery_saturation)

			aasRequested = np.floor(aasInSequences * translation_machinery_saturation)
		else:
			aasRequested = aasInSequences

		self.aas.requestIs(
			aasRequested
			)

		self.writeToListener("GrowthLimits", "aaPoolSize", self.aas.total())
		self.writeToListener("GrowthLimits", "aaRequestSize", aasRequested)

		if self.translationSaturation:
			gtpsHydrolyzed = np.int64(np.ceil(
				self.gtpPerElongation * np.fmin(
					sequenceHasAA.sum(),
					np.floor(self.aas.total().sum() * translation_machinery_saturation)
					)
				))
		else:
			gtpsHydrolyzed = np.int64(np.ceil(
				self.gtpPerElongation * np.fmin(
					sequenceHasAA.sum(),
					self.aas.total().sum()
					)
				))

		self.writeToListener("GrowthLimits", "gtpPoolSize", self.gtp.total()[0])
		self.writeToListener("GrowthLimits", "gtpRequestSize", gtpsHydrolyzed)

		self.gtp.requestIs(gtpsHydrolyzed)

		self.h2o.requestIs(gtpsHydrolyzed) # note: this is roughly a 2x overestimate


	# Calculate temporal evolution
	def evolveState(self):
		self.writeToListener("GrowthLimits", "gtpAllocated", self.gtp.count())
		self.writeToListener("GrowthLimits", "aaAllocated", self.aas.counts())

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

		# trnasCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.trna_groups],dtype = np.int64)
		# synthetaseCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.synthetase_groups],dtype = np.int64)
		# elongationResourceCapacity = np.minimum(aaCounts, synthetaseCapacity, trnasCapacity)

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
			(sequenceElongations > 0) &
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

		self.ribosome30S.countInc(nTerminated)
		self.ribosome50S.countInc(nTerminated)

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

		self.writeToListener("GrowthLimits", "aasUsed", aasUsed)
		self.writeToListener("GrowthLimits", "gtpUsed", gtpUsed)

		self.writeToListener("RibosomeData", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aaCounts)

		self.writeToListener("RibosomeData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RibosomeData", "actualElongations", sequenceElongations.sum())

		self.writeToListener("RibosomeData", "didTerminate", didTerminate.sum())
		self.writeToListener("RibosomeData", "terminationLoss", (terminalLengths - peptideLengths)[didTerminate].sum())
