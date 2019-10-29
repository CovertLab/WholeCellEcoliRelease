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
import copy
import re

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
from wholecell.utils.random import stochasticRound
from wholecell.utils import units



class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideElongation, self).initialize(sim, sim_data)

		# Load parameters
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity
		proteinIds = sim_data.process.translation.monomerData['id']
		self.proteinLengths = sim_data.process.translation.monomerData["length"].asNumber()
		self.proteinSequences = sim_data.process.translation.translationSequences
		self.aaWeightsIncorporated = sim_data.process.translation.translationMonomerWeights
		self.endWeight = sim_data.process.translation.translationEndWeight
		self.gtpPerElongation = sim_data.constants.gtpPerTranslation
		self.translation_data = sim_data.process.translation
		self.ribosomeElongationRate = float(sim_data.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s))

		self.base_elongation_rate = float(sim_data.constants.ribosomeElongationRateBase.asNumber(units.aa / units.s))

		self.ribosomeElongationRateDict = sim_data.process.translation.ribosomeElongationRateDict

		self.translation_aa_supply = sim_data.translationSupplyRate

		# Used for figure in publication
		self.trpAIndex = np.where(proteinIds == "TRYPSYN-APROTEIN[c]")[0][0]

		# Create view onto activly elongating 70S ribosomes
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		# Create views onto 30S and 50S ribosomal subunits for termination
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		# Create view onto all proteins
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		# Create views onto all polymerization reaction small molecules
		self.aas = self.bulkMoleculesView(sim_data.moleculeGroups.aaIDs)
		self.h2o = self.bulkMoleculeView('WATER[c]')
		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h   = self.bulkMoleculeView("PROTON[c]")

		# Set for timestep calculation
		self.gtpUsed = 0
		self.gtpAvailable = 0

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		self.translationSupply = sim._translationSupply
		self.flat_elongation = not sim._variable_elongation_translation

		# This gets set for daughters in initial_conditions.py
		self.elngRateFactor = 1.

	def calculateRequest(self):
		# Set ribosome elongation rate based on simulation medium environment and elongantion rate
		# factor which is used to create single-cell variability in growth rate. The maximum number
		# of amino acids that can be elongated in a single timestep is set to 22 intentionally as
		# the minimum number of padding values on the protein sequence matrix is set to 22. If
		# timesteps longer than 1.0s are used, this feature will lead to errors in the effective
		# ribosome elongation rate.

		current_nutrients = self._external_states['Environment'].nutrients

		if self.translationSupply:
			self.ribosomeElongationRate = self.base_elongation_rate
		else:
			rate = self.ribosomeElongationRateDict[current_nutrients].asNumber(units.aa / units.s)
			self.ribosomeElongationRate = self.elngRateFactor * rate

		self.elongation_rates = self.translation_data.make_elongation_rates(
			self.randomState,
			self.ribosomeElongationRate,
			self.timeStepSec(),
			self.flat_elongation)

		# Request all active ribosomes
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		if len(activeRibosomes) == 0:
			return

		# Build sequences to request appropriate amount of amino acids to
		# polymerize for next timestep
		proteinIndexes, peptideLengths = activeRibosomes.attrs(
			'proteinIndex',
			'peptideLength')

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elongation_rates)

		sequenceHasAA = (sequences != polymerize.PAD_VALUE)
		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		if self.translationSupply:
			translationSupplyRate = self.translation_aa_supply[current_nutrients] * self.elngRateFactor

			self.writeToListener("RibosomeData", "translationSupply", translationSupplyRate.asNumber())


			dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

			molAasRequested = translationSupplyRate * dryMass * self.timeStepSec() * units.s

			countAasRequested = units.convertNoUnitToNumber(molAasRequested * self.nAvogadro)

			countAasRequested = np.fmin(countAasRequested, aasInSequences) # Check if this is required. It is a better request but there may be fewer elongations.
		else:
			countAasRequested = aasInSequences

		self.aas.requestIs(
			countAasRequested
			)

		self.writeToListener("GrowthLimits", "aaPoolSize", self.aas.total())
		self.writeToListener("GrowthLimits", "aaRequestSize", countAasRequested)

		# Request GTP for polymerization based on sequences
		gtpsHydrolyzed = np.int64(np.ceil(self.gtpPerElongation * countAasRequested.sum()))

		self.writeToListener("GrowthLimits", "gtpPoolSize", self.gtp.total()[0])
		self.writeToListener("GrowthLimits", "gtpRequestSize", gtpsHydrolyzed)

		# GTP hydrolysis is carried out in Metabolism process for growth associated maintenence
		# THis is set here for metabolism to use
		self.gtpRequest = gtpsHydrolyzed

	def evolveState(self):
		# Write allocation data to listener
		self.writeToListener("GrowthLimits", "gtpAllocated", self.gtp.count())
		self.writeToListener("GrowthLimits", "aaAllocated", self.aas.counts())

		# Get number of active ribosomes
		activeRibosomes = self.activeRibosomes.molecules()

		self.writeToListener("GrowthLimits", "activeRibosomeAllocated", len(activeRibosomes))

		if len(activeRibosomes) == 0:
			return

		# Build amino acids sequences for each ribosome to polymerize
		proteinIndexes, peptideLengths, massDiffProtein = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength', 'massDiff_protein'
			)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elongation_rates)

		if sequences.size == 0:
			return

		# Calculate elongation resource capacity
		aaCountInSequence = np.bincount(sequences[(sequences != polymerize.PAD_VALUE)])
		aaCounts = self.aas.counts()

		# Using polymerization algorithm elongate each ribosome up to the limits
		# of amino acids, sequence, and GTP
		active_elongation_rates = self.elongation_rates[proteinIndexes]
		result = polymerize(
			sequences,
			aaCounts,
			10000000, # Set to a large number, the limit is now taken care of in metabolism
			self.randomState,
			active_elongation_rates)

		sequenceElongations = result.sequenceElongation
		aasUsed = result.monomerUsages
		nElongations = result.nReactions

		# Update masses of ribosomes attached to polymerizing polypeptides
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

		# Write current average elongation to listener
		currElongRate = (sequenceElongations.sum() / len(activeRibosomes)) / self.timeStepSec()
		self.writeToListener("RibosomeData", "effectiveElongationRate", currElongRate)

		# Update active ribosomes, terminating if neccessary
		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiff_protein = updatedMass
			)

		# Ribosomes that reach the end of their sequences are terminated and
		# dissociated into 30S and 50S subunits. The polypeptide that they are polymerizing
		# is converted into a protein in BulkMolecules
		terminalLengths = self.proteinLengths[proteinIndexes]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedProteins = np.bincount(
			proteinIndexes[didTerminate],
			minlength = self.proteinSequences.shape[0]
			)

		activeRibosomes.delByIndexes(np.where(didTerminate)[0])
		self.bulkMonomers.countsInc(terminatedProteins)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		self.ribosome30S.countInc(nTerminated)
		self.ribosome50S.countInc(nTerminated)

		# Update counts of amino acids and water to reflect polymerization reactions
		self.aas.countsDec(aasUsed)
		self.h2o.countInc(nElongations - nInitialized)

		# Report stalling information
		expectedElongations = np.fmin(
			self.ribosomeElongationRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		# Write data to listeners
		self.writeToListener("GrowthLimits", "aasUsed", aasUsed)
		self.writeToListener("GrowthLimits", "gtpUsed", self.gtpUsed)

		self.writeToListener("RibosomeData", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aaCounts)

		self.writeToListener("RibosomeData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RibosomeData", "actualElongations", sequenceElongations.sum())
		self.writeToListener("RibosomeData", "actualElongationHist", np.histogram(sequenceElongations, bins = np.arange(0,23))[0])
		self.writeToListener("RibosomeData", "elongationsNonTerminatingHist", np.histogram(sequenceElongations[~didTerminate], bins=np.arange(0,23))[0])

		self.writeToListener("RibosomeData", "didTerminate", didTerminate.sum())
		self.writeToListener("RibosomeData", "terminationLoss", (terminalLengths - peptideLengths)[didTerminate].sum())
		self.writeToListener("RibosomeData", "numTrpATerminated", terminatedProteins[self.trpAIndex])

		self.writeToListener("RibosomeData", "processElongationRate", self.ribosomeElongationRate / self.timeStepSec())

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		"""
		Assumes GTP is the readout for failed translation with respect to the timestep.
		"""

		# Until more padding values are added to the protein sequence matrix, limit the maximum timestep length to 1 second
		# Since the current upper limit on a.a's elongated by ribosomes during a single timestep is set to 22, timesteps
		# longer than 1.0s do not lead to errors, but does slow down the ribosome elongation rate of the resulting simulation.
		# Must be modified if timesteps longer than 1.0s are desired.
		if inputTimeStep > 1.0:
			return False

		activeRibosomes = float(self.activeRibosomes.total()[0])
		self.gtpAvailable = float(self.gtp.total()[0])

		# Without an estimate on ribosome counts, require a short timestep until estimates available
		if activeRibosomes == 0:
			if inputTimeStep <= .2:
				return True
			else:
				return False

		dt = inputTimeStep * timeStepSafetyFraction
		gtpExpectedUsage = activeRibosomes * self.ribosomeElongationRate * self.gtpPerElongation * dt

		if gtpExpectedUsage < self.gtpAvailable:
			return True
		else:
			return False

	def wasTimeStepShortEnough(self):
		"""
		If translation used more than 90 percent of gtp, timeStep was too short.
		"""

		# If gtpAvailable is 0 and the timeStep is short, use the gtp produced this timeStep as the estimate
		if self.gtpAvailable == 0 and self.timeStepSec() <= .2:
			self.gtpAvailable = self.gtp.total()[0]

		if (self.gtpAvailable * .9) < self.gtpUsed:
			return False
		else:
			return True
