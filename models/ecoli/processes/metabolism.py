#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- enzyme-limited reactions (& fit enzyme expression)
- option to call a reduced form of metabolism (assume optimal)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

from itertools import izip

import numpy as np
from scipy.sparse import csr_matrix

import wholecell.processes.process
from wholecell.utils import units

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

from wholecell.utils.modular_fba import FluxBalanceAnalysis
from wholecell.utils.enzymeKinetics import EnzymeKinetics

from wholecell.utils.fitting import massesAndCountsToAddForPools

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

SECRETION_PENALTY_COEFF = 1e-5

NONZERO_ENZYMES = False
USE_RATELIMITS = False # Enable/disable kinetic rate limits in the model
USE_BASE_RATES = False

USE_MANUAL_FLUX_COEFF = False # enable to overrid flux coefficients in the knowledgebase and use these local values instead
MAX_FLUX_COEFF = 2 # Multiple of predicted rate at which to set the max fluxes
MIN_FLUX_COEFF = 0 # Multiple of predicted rate at which to set the min fluxes


class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Metabolism, self).initialize(sim, sim_data)

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		self.metabolitePoolIDs = sorted(sim_data.process.metabolism.concDict)

		self.exchangeConstraints = sim_data.process.metabolism.exchangeConstraints

		self.doublingTime = sim_data.doubling_time
		self.nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel

		# Load enzyme kinetic rate information
		self.reactionRateInfo = sim_data.process.metabolism.reactionRateInfo
		self.enzymeNames = sim_data.process.metabolism.enzymeNames
		self.constraintIDs = sim_data.process.metabolism.constraintIDs
		self.constraintToReactionDict = sim_data.process.metabolism.constraintToReactionDict
		self.reactionEnzymes  = sim_data.process.metabolism.reactionEnzymes

		if USE_MANUAL_FLUX_COEFF:
			self.max_flux_coefficient = MAX_FLUX_COEFF
			self.min_flux_coefficient = MIN_FLUX_COEFF
		else:
			self.max_flux_coefficient = sim_data.constants.kineticRateLimitFactorUpper
			self.min_flux_coefficient = sim_data.constants.kineticRateLimitFactorLower

		self.objective = dict(
			(key, sim_data.process.metabolism.concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in sim_data.process.metabolism.concDict
			)

		self.getMass = sim_data.getter.getMass
		self.massReconstruction = sim_data.mass
		self.avgCellToInitialCellConvFactor = sim_data.mass.avgCellToInitialCellConvFactor

		self.ngam = sim_data.constants.nonGrowthAssociatedMaintenance

		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		self.energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		self.reactionStoich = sim_data.process.metabolism.reactionStoich
		self.externalExchangeMolecules = sim_data.nutrientData["secretionExchangeMolecules"]
		for time, nutrientsLabel in sim_data.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel]:
			self.externalExchangeMolecules += sim_data.nutrientData["importExchangeMolecules"][nutrientsLabel]
		self.externalExchangeMolecules = sorted(self.externalExchangeMolecules)
		self.extMoleculeMasses = self.getMass(self.externalExchangeMolecules)

		self.moleculeMasses = dict(zip(
			self.externalExchangeMolecules,
			self.getMass(self.externalExchangeMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)
			))

		# Set up FBA solver
		self.fba = FluxBalanceAnalysis(
			self.reactionStoich,
			self.externalExchangeMolecules,
			self.objective,
			objectiveType = "pools",
			moleculeMasses = self.moleculeMasses,
			secretionPenaltyCoeff = SECRETION_PENALTY_COEFF, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			solver = "glpk",
			maintenanceCostGAM = self.energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			maintenanceReaction = {
				"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "Pi[c]": +1, "PROTON[c]": +1,
				} # TODO: move to KB TODO: check reaction stoich
			)

		# Matrix mapping enzymes to the reactions they catalyze
		self.enzymeReactionMatrix = sim_data.process.metabolism.enzymeReactionMatrix(self.fba.reactionIDs(), self.enzymeNames, self.reactionEnzymes)
		self.spontaneousIndices = np.where(np.sum(self.enzymeReactionMatrix, axis=1) == 0)
		self.enzymeReactionMatrix = csr_matrix(self.enzymeReactionMatrix)
		self.kcat_max = sim_data.constants.carbonicAnhydraseKcat
		self.base_rates_current = FLUX_UNITS * np.inf * np.ones(len(self.fba.reactionIDs()))

		# # Set up enzyme kinetics object
		# self.enzymeKinetics = EnzymeKinetics(
		# 	enzymesWithKineticInfo = self.enzymeNames,
		# 	reactionRateInfo = self.reactionRateInfo,
		# 	constraintIDs = sim_data.process.metabolism.constraintIDs,
		# 	reactionIDs = self.fba.reactionIDs(),
		# 	metaboliteIDs = self.fba.outputMoleculeIDs(),
		# 	kcatOnly=False
		# 	)

		self.currentNgam = 1 * (COUNTS_UNITS / VOLUME_UNITS)
		self.currentPolypeptideElongationEnergy = 1 * (COUNTS_UNITS / VOLUME_UNITS)

		# # Determine which kinetic limits to use
		# self.reactionsWithKineticLimits = [True]*len(self.fba.reactionIDs())
	
		# Set constraints
		## External molecules
		self.externalMoleculeIDs = self.fba.externalMoleculeIDs()

		## Set enzymes unlimited
		self.fba.enzymeLevelsIs(np.inf)

		# Views
		self.metaboliteNames = self.fba.outputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNames)
		self.poolMetabolites = self.bulkMoleculesView(self.metabolitePoolIDs)
		self.enzymes = self.bulkMoleculesView(self.enzymeNames)
			
		outputMoleculeIDs = self.fba.outputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.internalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		self.fitterPredictedFluxesDict = sim_data.process.metabolism.predictedFluxesDict

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.enzymes.requestAll()

	def evolveState(self):

		# Solve for metabolic fluxes
		metaboliteCountsInit = self.metabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

		cellVolume = cellMass / self.cellDensity

		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		self.newPolypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			self.newPolypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest

		# Set external molecule levels
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)

		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			self.nutrientsTimeSeriesLabel,
			self.time()
			)

		if newObjective != None and newObjective != self.objective:
			# Build new fba instance with new objective
			self.objective = newObjective
			self.fba = FluxBalanceAnalysis(
				self.reactionStoich,
				self.externalExchangeMolecules,
				self.objective,
				objectiveType = "pools",
				moleculeMasses = self.moleculeMasses,
				secretionPenaltyCoeff = SECRETION_PENALTY_COEFF,
				solver = "glpk",
				maintenanceCostGAM = self.energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
				maintenanceReaction = {
					"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "Pi[c]": +1, "PROTON[c]": +1,
					} # TODO: move to KB TODO: check reaction stoich
				)

			massComposition = self.massReconstruction.getFractionMass(self.doublingTime)
			massInitial = (massComposition["proteinMass"] + massComposition["rnaMass"] + massComposition["dnaMass"]) / self.avgCellToInitialCellConvFactor
			objIds = sorted(self.objective)
			objConc = (COUNTS_UNITS / VOLUME_UNITS) * np.array([self.objective[x] for x in objIds])
			mws = self.getMass(objIds)
			massesToAdd, _ = massesAndCountsToAddForPools(massInitial, objIds, objConc, mws, self.cellDensity, self.nAvogadro)
			smallMoleculePoolsDryMass = units.hstack((massesToAdd[:objIds.index('WATER[c]')], massesToAdd[objIds.index('WATER[c]') + 1:]))
			totalDryMass = units.sum(smallMoleculePoolsDryMass) + massInitial
			self.writeToListener("CellDivision", "expectedDryMassIncrease", totalDryMass)

		# Set external molecule levels
		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		self.newNgam = self.ngam * coefficient

		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		ADJUSTMENT_RATIO = .01
		ngam_diff = np.abs(self.currentNgam.asNumber() - self.newNgam.asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		if ngam_diff > ADJUSTMENT_RATIO:
			self.currentNgam = self.newNgam
			self.fba.maxReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.fba.minReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))

		poly_diff = np.abs((self.currentPolypeptideElongationEnergy.asNumber() - self.newPolypeptideElongationEnergy.asNumber())) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		if poly_diff > ADJUSTMENT_RATIO:
			self.currentPolypeptideElongationEnergy = self.newPolypeptideElongationEnergy
			self.fba.maxReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.fba.minReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		#  Find metabolite concentrations from metabolite counts
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit

		self.fba.internalMoleculeLevelsIs(
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS)
			)

		#  Find enzyme concentrations from enzyme counts
		enzymeCountsInit = self.enzymes.counts()

		enzymeConcentrations = countsToMolar * enzymeCountsInit

		if NONZERO_ENZYMES:
			# Add one of every enzyme to ensure none are zero
			enzymeConcentrations = countsToMolar * (enzymeCountsInit + 1)

		if USE_BASE_RATES:
			# Calculate new rates
			self.enzymeMaxRates = self.kcat_max * enzymeConcentrations
			self.base_rates_new = FLUX_UNITS * self.enzymeReactionMatrix.dot(self.enzymeMaxRates.asNumber(FLUX_UNITS))
			self.base_rates_new[self.spontaneousIndices] = (FLUX_UNITS) * np.inf
			# Update base_rates_current
			updateLocations = np.where(self.base_rates_new.asNumber(FLUX_UNITS) != self.base_rates_current.asNumber(FLUX_UNITS))
			updateReactions = self.fba.reactionIDs()[updateLocations]
			updateValues = self.base_rates_new[updateLocations]
			self.base_rates_current[updateLocations] = updateValues
			# Set new reaction rate limits
			self.fba.setMaxReactionFluxes(updateReactions, updateValues.asNumber(FLUX_UNITS), raiseForReversible = False)

		# defaultRate = self.enzymeKinetics.defaultRate

		# # Combine the enzyme concentrations, substrate concentrations, and the default rate into one vector
		# inputConcentrations = np.concatenate((
		# 	enzymeConcentrations.asNumber(units.umol / units.L),
		# 	metaboliteConcentrations.asNumber(units.umol / units.L),
		# 	[defaultRate]), axis=1)

		# # Find reaction rate limits
		# self.reactionConstraints = ((units.umol / units.L) * self.enzymeKinetics.rateFunction(*inputConcentrations) * self.timeStepSec()).asNumber(COUNTS_UNITS / VOLUME_UNITS)

		# # Find rate limits for all constraints
		# self.allConstraintsLimits = ((units.umol / units.L) * self.enzymeKinetics.allRatesFunction(*inputConcentrations)[0]).asNumber(COUNTS_UNITS / VOLUME_UNITS)

		# # Set the rate limits only if the option flag is enabled
		# if USE_RATELIMITS:
		# 	currentRateLimits = {}
		# 	# Set reaction fluxes to be between  MAX_FLUX_COEFF and MIN_FLUX_COEFF of the predicted rate
		# 	for index, constraintID in enumerate(self.constraintIDs):
		# 		# Only use this kinetic limit if it's enabled
		# 		if self.activeConstraints[index]:
		# 			# Make sure to never set negative maximum rates
		# 			assert (self.allConstraintsLimits[index] >= 0 and self.allConstraintsLimits[index] != np.nan)

		# 			# Ensure that this reaction hasn't already been constrained more than this yet
		# 			if self.constraintToReactionDict[constraintID] in currentRateLimits and currentRateLimits[self.constraintToReactionDict[constraintID]] < self.allConstraintsLimits[index]*MAX_FLUX_COEFF:
		# 				# This rate has already been constrained more than this constraint, so skip it
		# 				continue

		# 			# Set the max reaction rate for this reaction
		# 			self.fba.maxReactionFluxIs(self.constraintToReactionDict[constraintID], self.allConstraintsLimits[index]*self.max_flux_coefficient, raiseForReversible = False)
		# 			# Set the minimum reaction rate for this reaction
		# 			self.fba.minReactionFluxIs(self.constraintToReactionDict[constraintID], self.allConstraintsLimits[index]*self.min_flux_coefficient, raiseForReversible = False)
					
		# 			# Record what constraint was just applied to this reaction
		# 			currentRateLimits[self.constraintToReactionDict[constraintID]] = self.allConstraintsLimits[index]*MAX_FLUX_COEFF

		# 		else:
		# 			self.fba.maxReactionFluxIs(self.constraintToReactionDict[constraintID], defaultRate, raiseForReversible = False)

		# self.overconstraintMultiples = (self.fba.reactionFluxes() / self.timeStepSec()) / self.reactionConstraints

		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.outputMoleculeLevelsChange())

		metaboliteCountsFinal = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.externalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)

		# TODO: report as reactions (#) per second & store volume elsewhere
		self.writeToListener("FBAResults", "reactionFluxes",
			self.fba.reactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes",
			exFluxes)
		# self.writeToListener("FBAResults", "objectiveValue", # TODO
		# 	self.fba.objectiveValue() / deltaMetabolites.size) # divide to normalize by number of metabolites
		self.writeToListener("FBAResults", "outputFluxes",
			self.fba.outputMoleculeLevelsChange() / self.timeStepSec())

		self.writeToListener("FBAResults", "rowDualValues",
			self.fba.rowDualValues(self.metaboliteNames))

		self.writeToListener("FBAResults", "columnDualValues",
			self.fba.columnDualValues(self.fba.reactionIDs()))

		# self.writeToListener("EnzymeKinetics", "reactionConstraints",
		# 	self.reactionConstraints)

		# self.writeToListener("EnzymeKinetics", "allConstraintsLimits",
		# 	self.allConstraintsLimits)

		# self.writeToListener("EnzymeKinetics", "overconstraintMultiples",
		# 	self.overconstraintMultiples)

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit",
			metaboliteCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal",
			metaboliteCountsFinal)

		self.writeToListener("EnzymeKinetics", "enzymeCountsInit",
			enzymeCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations",
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		self.writeToListener("EnzymeKinetics", "countsToMolar",
			countsToMolar.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		self.writeToListener("EnzymeKinetics", "counts_units",
			str(COUNTS_UNITS))

		self.writeToListener("EnzymeKinetics", "mass_units",
			str(MASS_UNITS))

		self.writeToListener("EnzymeKinetics", "volume_units",
			str(VOLUME_UNITS))




		# TODO
		# NOTE: the calculation for the objective components doesn't yet have
		# an interface, since it will vary in calculation and shape for every
		# objective type

		# objectiveComponents_raw = (np.array(self.fba._f).flatten() * self.fba._solutionFluxes)[self.fba._objIndexes]
		# objectiveComponents = objectiveComponents_raw[::2] + objectiveComponents_raw[1::2]

		# self.writeToListener("FBAResults", "objectiveComponents",
		# 	objectiveComponents
		# 	)

		# TODO:
		# - which media exchanges/reactions are limiting, if any
		# - objective details (value, component values)
