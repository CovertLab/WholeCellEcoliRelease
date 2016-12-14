#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
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

from wholecell.utils.fitting import massesAndCountsToAddForHomeostaticTargets

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

SECRETION_PENALTY_COEFF = 1e-5

NONZERO_ENZYMES = False

USE_KINETIC_RATES = True
USE_BASE_RATES = True
KINETICS_BURN_IN_PERIOD = 1

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

		self.exchangeConstraints = sim_data.process.metabolism.exchangeConstraints

		self.doublingTime = sim_data.doubling_time
		self.nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel

		# Load enzyme kinetic rate information
		self.reactionRateInfo = sim_data.process.metabolism.reactionRateInfo
		self.enzymeNames = sim_data.process.metabolism.metabolicEnzymes
		self.constraintIDs = sim_data.process.metabolism.constraintIDs
		self.constraintMultiplesDict = {constraintID:rateInfo["constraintMultiple"] for constraintID, rateInfo in self.reactionRateInfo.iteritems()}
		self.constraintToReactionDict = sim_data.process.metabolism.constraintToReactionDict
		self.reactionEnzymes  = sim_data.process.metabolism.reactionEnzymes
		self.getBiomassAsConcentrations = sim_data.mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel][0][1]
			)

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[sim_data.condition])
		concDict.update(self.concModificationsBasedOnCondition)

		self.objective = dict(
			(key, concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in concDict
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

		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel]:
			self.externalExchangeMolecules += sim_data.nutrientData["importExchangeMolecules"][nutrientsLabel]

			# Sorry
			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					nutrientsLabel, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)

		self.maintenanceReaction = sim_data.process.metabolism.maintenanceReaction
		self.externalExchangeMolecules = sorted(self.externalExchangeMolecules)
		self.extMoleculeMasses = self.getMass(self.externalExchangeMolecules)

		self.moleculeMasses = dict(zip(
			self.externalExchangeMolecules,
			self.getMass(self.externalExchangeMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)
			))

		# Set up enzyme kinetics object
		self.enzymeKinetics = EnzymeKinetics(
			reactionRateInfo = sim_data.process.metabolism.reactionRateInfo,
			useCustoms=True,
			moreThanKcat=False, # Only compute rates for reactions with more than a kcat
		)

		# Remove kinetics for reactions for which we don't have needed metabolites or enzymes
		self.enzymeKinetics.checkKnownSubstratesAndEnzymes(sim_data.process.metabolism.concDict, self.enzymeNames, removeUnknowns=True)

		# Add reactions with a kinetic estimate
		self.allRateReactions = sorted(set([reactionInfo["reactionID"] for constraintID, reactionInfo in self.enzymeKinetics.reactionRateInfo.iteritems() if reactionInfo["reactionID"] in self.reactionStoich]))
		# Reactions with full kinetic estimates (more than just kcat)
		self.fullRateReactions = sorted(set([reactionInfo["reactionID"] for constraintID, reactionInfo in self.enzymeKinetics.reactionRateInfo.iteritems() if (len(reactionInfo["kM"]) > 0 or reactionInfo["rateEquationType"] == "custom") and reactionInfo["reactionID"] in self.reactionStoich]))
		# Reactions with a kcat-based kinetic estimate only (no customs, no kMs, no kIs)
		self.kcatRateReactions = sorted(set([reactionInfo["reactionID"] for constraintID, reactionInfo in self.enzymeKinetics.reactionRateInfo.iteritems() if reactionInfo["reactionID"] not in self.fullRateReactions and reactionInfo["reactionID"] in self.reactionStoich]))

		self.metabolismKineticObjectiveWeight = sim_data.constants.metabolismKineticObjectiveWeight

		# Set up FBA solver
		self.fbaObjectOptions = {
			"reactionStoich" : self.reactionStoich,
			"externalExchangedMolecules" : self.externalExchangeMolecules,
			"objective" : self.objective,
			"objectiveType" : "homeostatic_kinetics_mixed",
			"objectiveParameters" : {
					"kineticObjectiveWeight":self.metabolismKineticObjectiveWeight,
					"reactionRateTargets":{reaction:1e-5 for reaction in self.allRateReactions}, #This target is arbitrary, it gets reset each timestep during evolveState
					"oneSidedReactionTargets":self.kcatRateReactions,
					},
			"moleculeMasses" : self.moleculeMasses,
			"secretionPenaltyCoeff" : SECRETION_PENALTY_COEFF, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver" : "glpk",
			"maintenanceCostGAM" : self.energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction" : self.maintenanceReaction,
		}

		if USE_KINETIC_RATES==False:
			self.fbaObjectOptions["objectiveType"] = "homeostatic"

		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.outputMoleculeIDs()])

		# Disable all rates during burn-in
		if KINETICS_BURN_IN_PERIOD > 0 and USE_KINETIC_RATES:
			self.fba.disableKineticTargets()
			self.burnInComplete = False
		else:
			self.burnInComplete = True

		# Indices for reactions with full kinetic esimates
		self.allRateIndices = np.where([True if reactionID in self.allRateReactions else False for reactionID in self.fba.reactionIDs()])

		self.allRateEstimates = FLUX_UNITS * np.zeros_like(self.allRateIndices)

		# Matrix mapping enzymes to the reactions they catalyze
		self.enzymeReactionMatrix = sim_data.process.metabolism.enzymeReactionMatrix(
			self.fba.reactionIDs(),
			self.enzymeNames,
			self.reactionEnzymes,
			)
		self.spontaneousIndices = np.where(np.sum(self.enzymeReactionMatrix, axis=1) == 0)
		self.enzymeReactionMatrix = csr_matrix(self.enzymeReactionMatrix)
		self.baseRates = FLUX_UNITS * np.inf * np.ones(len(self.fba.reactionIDs()))

		self.currentNgam = 1 * (COUNTS_UNITS / VOLUME_UNITS)
		self.currentPolypeptideElongationEnergy = 1 * (COUNTS_UNITS / VOLUME_UNITS)

		# Set constraints
		## External molecules
		self.externalMoleculeIDs = self.fba.externalMoleculeIDs()

		# Views
		self.metaboliteNames = self.fba.outputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
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

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(self.nutrientToDoublingTime.get(self._sim.processes["PolypeptideElongation"].currentNutrients, self.nutrientToDoublingTime["minimal"]))

		# Set external molecule levels
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)

		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			self.nutrientsTimeSeriesLabel,
			self.time(),
			self.concModificationsBasedOnCondition,
			)

		if newObjective != None and newObjective != self.objective:
			# Build new fba instance with new objective
			self.fbaObjectOptions["objective"] = newObjective
			self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
			self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.outputMoleculeIDs()])

		# After completing the burn-in, enable kinetic rates
		if self._sim.time() - self._sim.initialTime() > KINETICS_BURN_IN_PERIOD and USE_KINETIC_RATES and not self.burnInComplete:
			self.burnInComplete = True
			self.fba.enableKineticTargets()

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
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit[self.internalExchangeIdxs]

		# Make a dictionary of metabolite names to metabolite concentrations
		metaboliteConcentrationsDict = dict(zip(self.metaboliteNames, metaboliteConcentrations))

		self.fba.internalMoleculeLevelsIs(
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS)
			)

		#  Find enzyme concentrations from enzyme counts
		enzymeCountsInit = self.enzymes.counts()

		enzymeConcentrations = countsToMolar * enzymeCountsInit

		if NONZERO_ENZYMES:
			# Add one of every enzyme to ensure none are zero
			enzymeConcentrations = countsToMolar * (enzymeCountsInit + 1)

		# Make a dictionary of enzyme names to enzyme concentrations
		enzymeConcentrationsDict = dict(zip(self.enzymeNames, enzymeConcentrations))

		# When many estimates exist for a reaction, choose the largest
		if not hasattr(self, "chosenConstraints"):
			# Calculate the constraints in the current conditions
			reactionsDict = self.enzymeKinetics.allReactionsDict(metaboliteConcentrationsDict, enzymeConcentrationsDict)
			oneSidedReactions =  set(self.fba.kineticOneSidedTargetFluxNames())
			self.chosenConstraints = {}
			for reactionID, reactionRate in reactionsDict.iteritems():
				rateOrderedConstraints = sorted(reactionRate.keys(), key=reactionRate.__getitem__, reverse=False)
				kMreactions = [x for x in rateOrderedConstraints if 'kcat' not in x]
				if len(kMreactions) > 0:
					# Take the highest valued constraint with a kM
					constraintID = kMreactions[-1]
				elif len(kMreactions) == 0:
					# Take the higest valued constraint overall
					constraintID = rateOrderedConstraints[-1]
				self.chosenConstraints[reactionID] = {
					"constraintID":constraintID,
					"coefficient":self.constraintMultiplesDict[constraintID],}

		if USE_KINETIC_RATES and self.burnInComplete:
			self.allRateEstimates = self.enzymeKinetics.ratesView(self.allRateReactions, self.chosenConstraints, metaboliteConcentrationsDict, enzymeConcentrationsDict, raiseIfNotFound=True)

			# Make kinetic targets numerical zero instead of actually zero for solver stability
			self.allRateEstimates[self.allRateEstimates.asNumber() == 0] = FLUX_UNITS * 1e-20
			self.fba.setKineticTarget(self.allRateReactions, (TIME_UNITS*self.timeStepSec()*self.allRateEstimates).asNumber(COUNTS_UNITS/VOLUME_UNITS), raiseForReversible=False)

		if USE_BASE_RATES:
			# Calculate new rates
			self.baseRatesNew = FLUX_UNITS * self.enzymeReactionMatrix.dot(enzymeConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.baseRatesNew[self.spontaneousIndices] = (FLUX_UNITS) * np.inf
			# Update allRates
			updateLocations = np.where(self.baseRatesNew.asNumber(FLUX_UNITS) != self.baseRates.asNumber(FLUX_UNITS))
			updateReactions = self.fba.reactionIDs()[updateLocations]
			updateValues = self.baseRatesNew[updateLocations]
			self.baseRates[updateLocations] = updateValues
			# Set new reaction rate limits
			self.fba.setMaxReactionFluxes(updateReactions, (TIME_UNITS*self.timeStepSec()*updateValues).asNumber(COUNTS_UNITS/VOLUME_UNITS), raiseForReversible = False)

		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.outputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		# Use GLPK's dualprimal solver, AFTER the first solution
		self.fba._solver._model.set_solver_method_dualprimal()

		self.overconstraintMultiples = (self.fba.reactionFluxes()[self.allRateIndices] / self.allRateEstimates.asNumber(FLUX_UNITS))  
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

		self.writeToListener("FBAResults", "kineticObjectiveValues",
			self.fba.kineticObjectiveValues())

		self.writeToListener("FBAResults", "homeostaticObjectiveValues",
			self.fba.homeostaticObjectiveValues())

		self.writeToListener("FBAResults", "homeostaticObjectiveWeight",
			self.fba.homeostaticObjectiveWeight())

		self.writeToListener("EnzymeKinetics", "baseRates",
			self.baseRates.asNumber(FLUX_UNITS))

		self.writeToListener("EnzymeKinetics", "reactionKineticPredictions",
			self.allRateEstimates.asNumber(FLUX_UNITS))

		self.writeToListener("EnzymeKinetics", "overconstraintMultiples",
			self.overconstraintMultiples)

		self.writeToListener("EnzymeKinetics", "kineticTargetFluxes",
			self.fba.kineticTargetFluxes())

		self.writeToListener("EnzymeKinetics", "kineticTargetErrors",
			self.fba.kineticTargetFluxErrors())

		self.writeToListener("EnzymeKinetics", "kineticTargetRelativeDifferences",
			self.fba.kineticTargetFluxRelativeDifferences())

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
