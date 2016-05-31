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

# SYNTHETASE_KM_SCALE = 0.7

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		# Parameters
		self.maxRibosomeElongationRate = None
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
	def initialize(self, sim, sim_data):
		super(PolypeptideElongation, self).initialize(sim, sim_data)

		# Load parameters
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		proteinIds = sim_data.process.translation.monomerData['id']

		self.proteinLengths = sim_data.process.translation.monomerData["length"].asNumber()

		self.proteinSequences = sim_data.process.translation.translationSequences

		self.aaWeightsIncorporated = sim_data.process.translation.translationMonomerWeights

		self.endWeight = sim_data.process.translation.translationEndWeight

		self.gtpPerElongation = sim_data.constants.gtpPerTranslation

		self.maxRibosomeElongationRate = float(sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

		##########
		aaConcentrations = units.getUnit(sim_data.process.metabolism.concDict.items()[0][1]) * np.array([sim_data.process.metabolism.concDict[x].asNumber() for x in sim_data.moleculeGroups.aaIDs])
		total_aa_concentration = units.sum(aaConcentrations)
		sim_data.synthetase_km_scale = 0.3
		self.saturation_km = sim_data.synthetase_km_scale * total_aa_concentration
		##########

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		self.aas = self.bulkMoleculesView(sim_data.moleculeGroups.aaIDs)
		self.h2o = self.bulkMoleculeView('WATER[c]')

		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("Pi[c]")
		self.h   = self.bulkMoleculeView("PROTON[c]")

		self.gtpUsed = 0
		self.gtpAvailable = 0

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeGroups.s30_fullComplex[0])
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeGroups.s50_fullComplex[0])

		###### VARIANT CODE #######
		self.translationSaturation = sim_data.translationSaturation
		self.synthetase_km_scale = sim_data.synthetase_km_scale
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
			self._elngRate()
			)

		sequenceHasAA = (sequences != PAD_VALUE)

		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		if self.translationSaturation:
			cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
			cellVolume = cellMass / self.cellDensity
			aaTotalConc = units.sum((1 / self.nAvogadro) * (1 / cellVolume) * self.aas.total())

			# translation_machinery_saturation = (1 + self.synthetase_km_scale) * (aaTotalConc / (self.saturation_km + aaTotalConc))
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

		self.gtpRequest = gtpsHydrolyzed
		# self.gtp.requestIs(gtpsHydrolyzed)

		# self.h2o.requestIs(gtpsHydrolyzed) # note: this is roughly a 2x overestimate

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
			self._elngRate()
			)

		# Calculate elongation resource capacity
		aaCountInSequence = np.bincount(sequences[(sequences != PAD_VALUE)])
		aaCounts = self.aas.counts()

		# Calculate update
		reactionLimit = self.gtp.count() // self.gtpPerElongation

		sequenceElongations, aasUsed, nElongations = polymerize(
			sequences,
			aaCounts, # elongationResourceCapacity,
			10000000,#reactionLimit,
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

		# Write current average elongation rate for growth rate control
		currElongRate = (sequenceElongations.sum() / len(activeRibosomes)) / self.timeStepSec()
		self.writeToListener("RibosomeData", "effectiveElongationRate", currElongRate)
		print "Elng rate: {}".format(currElongRate)
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

		self.gtpUsed = 0#np.int64(stochasticRound(
		# 	self.randomState,
		# 	nElongations * self.gtpPerElongation
		# 	))

		self.gtp.countDec(self.gtpUsed)
		self.gdp.countInc(self.gtpUsed)
		self.pi.countInc(self.gtpUsed)
		self.h.countInc(self.gtpUsed)

		self.h2o.countDec(self.gtpUsed)

		# Report stalling information

		expectedElongations = np.fmin(
			self._elngRate(),
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		self.writeToListener("GrowthLimits", "aasUsed", aasUsed)
		self.writeToListener("GrowthLimits", "gtpUsed", self.gtpUsed)

		self.writeToListener("RibosomeData", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aaCounts)

		self.writeToListener("RibosomeData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RibosomeData", "actualElongations", sequenceElongations.sum())

		self.writeToListener("RibosomeData", "didTerminate", didTerminate.sum())
		self.writeToListener("RibosomeData", "terminationLoss", (terminalLengths - peptideLengths)[didTerminate].sum())

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		"""
		Assumes GTP is the readout for failed translation with respect to the timestep.
		"""

		activeRibosomes = float(self.activeRibosomes.total()[0])
		self.gtpAvailable = float(self.gtp.total()[0])

		# Without an estimate on ribosome counts, require a short timestep until estimates available
		if activeRibosomes == 0:
			if inputTimeStep <= .2:
				return True
			else:
				return False

		dt = inputTimeStep * timeStepSafetyFraction
		gtpExpectedUsage = activeRibosomes * self.maxRibosomeElongationRate * self.gtpPerElongation * dt

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

	def _elngRate(self):
		return int(stochasticRound(self.randomState, self.maxRibosomeElongationRate * self.timeStepSec()))
