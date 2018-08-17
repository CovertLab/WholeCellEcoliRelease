#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- option to call a reduced form of metabolism (assume optimal)
- handle oneSidedReaction constraints

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np
from scipy.sparse import csr_matrix

import wholecell.processes.process
from wholecell.utils import units

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

from wholecell.utils.modular_fba import FluxBalanceAnalysis

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

USE_KINETICS = True
KINETICS_BURN_IN_PERIOD = 0

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
		self.ngam = sim_data.constants.nonGrowthAssociatedMaintenance

		self.exchangeConstraints = sim_data.process.metabolism.exchangeConstraints
		self.updateExchangeData = sim_data.process.metabolism.exchangeDataFromConcentrations

		self.getBiomassAsConcentrations = sim_data.mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		# Create objective for homeostatic constraints
		nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label
		initial_environment = sim_data.external_state.environment.nutrients_time_series[nutrients_time_series_label][0][1]



		# initialize exchange_data according to initial concentrations in environment
		self.exchange_data = self.updateExchangeData(sim_data.external_state.environment.environment_dict[initial_environment])

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			initial_environment
			)
		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			sim_data.conditionToDoublingTime[sim_data.condition]
			)
		concDict.update(self.concModificationsBasedOnCondition)
		self.homeostaticObjective = dict((key, concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in concDict)

		# Load initial mass
		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit
		initCellMass = initWaterMass + initDryMass

		energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		# Setup molecules in external environment that can be exchanged
		#TODO (Eran) this can be replaced with reference to exchange_data
		externalExchangedMolecules = sim_data.process.metabolism.exchange_data_dict["secretionExchangeMolecules"][:]
		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.external_state.environment.nutrients_time_series[nutrients_time_series_label]:
			externalExchangedMolecules += sim_data.process.metabolism.exchange_data_dict["importExchangeMolecules"][nutrientsLabel]

			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					nutrientsLabel, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		externalExchangedMolecules = sorted(set(externalExchangedMolecules))

		# save nutrient names for environment view, using all moleculeIDs in local environment
		# TODO (Eran) make all environment molecules into external exchange molecules, not just this subset
		self.environment_molecule_ids = self._external_states['Environment']._moleculeIDs
		self.external_exchange_molecule_ids = externalExchangedMolecules

		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)

		moleculeMasses = dict(zip(externalExchangedMolecules, sim_data.getter.getMass(externalExchangedMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)))

		# Data structures to compute reaction bounds based on enzyme presence/absence
		self.catalystsList = sim_data.process.metabolism.catalystsList
		self.reactionsWithCatalystsList = sim_data.process.metabolism.reactionCatalystsList
		self.reactionCatalystsDict = sim_data.process.metabolism.reactionCatalysts

		catalysisMatrixI = sim_data.process.metabolism.catalysisMatrixI
		catalysisMatrixJ = sim_data.process.metabolism.catalysisMatrixJ
		catalysisMatrixV = sim_data.process.metabolism.catalysisMatrixV

		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
		self.catalysisMatrix = csr_matrix((catalysisMatrixV, (catalysisMatrixI, catalysisMatrixJ)), shape = shape)

		self.catalyzedReactionBoundsPrev = np.inf * np.ones(len(self.reactionsWithCatalystsList))

		# Function to compute reaction targets based on kinetic parameters and molecule concentrations
		self.getKineticConstraints = sim_data.process.metabolism.getKineticConstraints

		# Remove disabled reactions so they don't get included in the FBA problem setup
		if hasattr(sim_data.process.metabolism, "kineticTargetShuffleRxns") and sim_data.process.metabolism.kineticTargetShuffleRxns != None:
			self.kineticsConstrainedReactions = sim_data.process.metabolism.kineticTargetShuffleRxns
			self.active_constraints_mask = np.ones(len(self.kineticsConstrainedReactions), dtype=bool)
		else:
			constrainedReactionList = sim_data.process.metabolism.constrainedReactionList
			constraintsToDisable = sim_data.process.metabolism.constraintsToDisable
			self.active_constraints_mask = np.array([(rxn not in constraintsToDisable) for rxn in constrainedReactionList])
			self.kineticsConstrainedReactions = list(np.array(constrainedReactionList)[self.active_constraints_mask])

		self.kineticsEnzymesList = sim_data.process.metabolism.enzymeIdList
		self.kineticsSubstratesList = sim_data.process.metabolism.kineticsSubstratesList

		constraintToReactionMatrixI = sim_data.process.metabolism.constraintToReactionMatrixI
		constraintToReactionMatrixJ = sim_data.process.metabolism.constraintToReactionMatrixJ
		constraintToReactionMatrixV = sim_data.process.metabolism.constraintToReactionMatrixV
		shape = (constraintToReactionMatrixI.max() + 1, constraintToReactionMatrixJ.max() + 1)
		self.constraintToReactionMatrix = np.zeros(shape, np.float64)
		self.constraintToReactionMatrix[constraintToReactionMatrixI, constraintToReactionMatrixJ] = constraintToReactionMatrixV
		self.constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly

		# Set solver and kinetic objective weight (lambda)
		solver = sim_data.process.metabolism.solver
		kinetic_objective_weight = sim_data.process.metabolism.kinetic_objective_weight

		# Disable kinetics completely if weight is 0 or specified in file above
		self.use_kinetics = True
		if not USE_KINETICS or kinetic_objective_weight == 0:
			self.use_kinetics = False
			kinetic_objective_weight = 0

		# Set up FBA solver
		# reactionRateTargets value is just for initialization, it gets reset each timestep during evolveState
		self.fbaObjectOptions = {
			"reactionStoich" : sim_data.process.metabolism.reactionStoich,
			"externalExchangedMolecules" : externalExchangedMolecules,
			"objective" : self.homeostaticObjective,
			"objectiveType" : "homeostatic_kinetics_mixed",
			"objectiveParameters" : {
					"kineticObjectiveWeight" : kinetic_objective_weight,
					"reactionRateTargets" : {reaction : 1 for reaction in self.kineticsConstrainedReactions},
					"oneSidedReactionTargets" : [],
					},
			"moleculeMasses" : moleculeMasses,
			"secretionPenaltyCoeff" : sim_data.constants.secretion_penalty_coeff, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver" : solver,
			"maintenanceCostGAM" : energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction" : sim_data.process.metabolism.maintenanceReaction,
		}
		if not self.use_kinetics:
			self.fbaObjectOptions["objectiveType"] = "homeostatic"
		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])

		# Disable all rates during burn-in
		if self.use_kinetics:
			if KINETICS_BURN_IN_PERIOD > 0:
				self.fba.disableKineticTargets()
				self.burnInComplete = False
			else:
				self.burnInComplete = True

		# Values will get updated at each time point
		self.currentNgam = 1 * (COUNTS_UNITS / VOLUME_UNITS)
		self.currentPolypeptideElongationEnergy = 1 * (COUNTS_UNITS / VOLUME_UNITS)

		# External molecules
		self.externalMoleculeIDs = self.fba.getExternalMoleculeIDs()

		## Construct views
		# views on environment
		self.environment_molecules = self.environmentView(self.environment_molecule_ids)
		self.external_exchange_molecules = self.environmentView(self.external_exchange_molecule_ids)

		# views on metabolism bulk molecules
		self.metaboliteNames = self.fba.getOutputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.catalystsList)
		self.kineticsEnzymes = self.bulkMoleculesView(self.kineticsEnzymesList)
		self.kineticsSubstrates = self.bulkMoleculesView(self.kineticsSubstratesList)

		outputMoleculeIDs = self.fba.getOutputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.getInternalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		self.AAs = [x[:-3] for x in sorted(sim_data.amino_acid_1_to_3_ordered.values())]

		self.shuffleIdxs = None
		if hasattr(sim_data.process.metabolism, "kineticTargetShuffleIdxs") and sim_data.process.metabolism.kineticTargetShuffleIdxs != None:
			self.shuffleIdxs = sim_data.process.metabolism.kineticTargetShuffleIdxs

		self.shuffleCatalyzedIdxs = None
		if hasattr(sim_data.process.metabolism, "catalystShuffleIdxs") and sim_data.process.metabolism.catalystShuffleIdxs != None:
			self.shuffleCatalyzedIdxs = sim_data.process.metabolism.catalystShuffleIdxs

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.catalysts.requestAll()
		self.kineticsEnzymes.requestAll()
		self.kineticsSubstrates.requestAll()

	def evolveState(self):
		metaboliteCountsInit = self.metabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

		cellVolume = cellMass / self.cellDensity
		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		# TODO (Eran) remove dependence on current_nutrient label, work only with current_environment's concentrations
		current_nutrients = self._external_states['Environment'].nutrients

		# recalculate exchange_data based on current environment
		current_environment = dict(zip(self.environment_molecule_ids, self.environment_molecules.totalConcentrations() * COUNTS_UNITS / VOLUME_UNITS))
		self.exchange_data = self.updateExchangeData(current_environment)

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			self.nutrientToDoublingTime.get(current_nutrients, self.nutrientToDoublingTime["minimal"])
			)

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)

		# Set external molecule levels
		# TODO (Eran) remove current_nutrients
		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			current_nutrients,
			self.exchange_data,
			self.concModificationsBasedOnCondition,
			)

		updatedObjective = False
		if newObjective != None and newObjective != self.homeostaticObjective:
			# Build new fba instance with new objective
			self.fbaObjectOptions["objective"] = newObjective
			self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
			self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])
			self.homeostaticObjective = newObjective
			updatedObjective = True

		# After completing the burn-in, enable kinetic rates
		if self.use_kinetics and (not self.burnInComplete) and (self._sim.time() > KINETICS_BURN_IN_PERIOD):
			self.burnInComplete = True
			self.fba.enableKineticTargets()

		# Allow flexibility for solver in first time step after an environment shift
		if updatedObjective:
			self.fba.disableKineticTargets()
			self.burnInComplete = False

		#  Find metabolite concentrations from metabolite counts
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit[self.internalExchangeIdxs]

		# Make a dictionary of metabolite names to metabolite concentrations
		metaboliteConcentrationsDict = dict(zip(self.metaboliteNames, metaboliteConcentrations))
		self.fba.setInternalMoleculeLevels(metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# Set external molecule levels
		self._setExternalMoleculeLevels(externalMoleculeLevels, metaboliteConcentrations)

		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		ADJUSTMENT_RATIO = .01

		# Calculate new NGAM and update if necessary
		self.newNgam = self.ngam * coefficient
		ngam_diff = np.abs(self.currentNgam.asNumber() - self.newNgam.asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		if ngam_diff > ADJUSTMENT_RATIO:
			self.currentNgam = self.newNgam
			flux = (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS)
			self.fba.setReactionFluxBounds(self.fba._reactionID_NGAM, lowerBounds=flux, upperBounds=flux)

		# Calculate GTP usage based on how much was needed in polypeptide elongation in previous step and update if necessary
		newPolypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			newPolypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest
		poly_diff = np.abs((self.currentPolypeptideElongationEnergy.asNumber() - newPolypeptideElongationEnergy.asNumber())) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		if poly_diff > ADJUSTMENT_RATIO:
			self.currentPolypeptideElongationEnergy = newPolypeptideElongationEnergy
			flux = self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS)
			self.fba.setReactionFluxBounds(self.fba._reactionID_polypeptideElongationEnergy, lowerBounds=flux, upperBounds=flux)

		# Constrain reactions based on absence of catalysts
		## Read counts for catalysts and enzymes (catalysts with kinetics constraints)
		catalystsCountsInit = self.catalysts.counts()

		## Set hard upper bounds constraints based on enzyme presence (infinite upper bound) or absence (upper bound of zero)
		catalyzedReactionBounds = np.inf * np.ones(len(self.reactionsWithCatalystsList))
		rxnPresence = self.catalysisMatrix.dot(catalystsCountsInit)
		catalyzedReactionBounds[rxnPresence == 0] = 0
		if self.shuffleCatalyzedIdxs is not None:
			catalyzedReactionBounds = catalyzedReactionBounds[self.shuffleCatalyzedIdxs]

		## Only update reaction limits that are different from previous time step
		updateIdxs = np.where(catalyzedReactionBounds != self.catalyzedReactionBoundsPrev)[0]
		updateRxns = [self.reactionsWithCatalystsList[idx] for idx in updateIdxs]
		updateVals = catalyzedReactionBounds[updateIdxs]
		self.fba.setReactionFluxBounds(updateRxns, upperBounds=updateVals, raiseForReversible=False)
		self.catalyzedReactionBoundsPrev = catalyzedReactionBounds

		# Constrain reactions based on kinetic values
		kineticsEnzymesCountsInit = self.kineticsEnzymes.counts()
		kineticsEnzymesConcentrations = countsToMolar * kineticsEnzymesCountsInit

		kineticsSubstratesCountsInit = self.kineticsSubstrates.counts()
		kineticsSubstratesConcentrations = countsToMolar * kineticsSubstratesCountsInit

		## Set target fluxes for reactions based on their most relaxed constraint
		constraintValues = self.getKineticConstraints(
			kineticsEnzymesConcentrations.asNumber(units.umol / units.L),
			kineticsSubstratesConcentrations.asNumber(units.umol / units.L),
			)
		reactionTargets = (units.umol / units.L / units.s) * np.max(self.constraintToReactionMatrix * constraintValues, axis = 1)

		## Shuffle parameters (only performed in very specific cases)
		if self.shuffleIdxs is not None:
			reactionTargets = (units.umol / units.L / units.s) * reactionTargets.asNumber()[self.shuffleIdxs]

		## Record which constraint was used, add constraintToReactionMatrix to ensure the index is one of the constraints if multiplication is 0
		reactionConstraint = np.argmax(self.constraintToReactionMatrix * constraintValues + self.constraintToReactionMatrix, axis = 1)

		## Calculate reaction flux target for current time step
		targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(COUNTS_UNITS / VOLUME_UNITS)[self.active_constraints_mask]

		## Set kinetic targets only if kinetics is enabled
		if self.use_kinetics and self.burnInComplete:
			self.fba.setKineticTarget(self.kineticsConstrainedReactions, targets, raiseForReversible = False)

		# Solve FBA problem and update metabolite counts
		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.getOutputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.getExternalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)

		# change in nutrient counts, used in non-infinite environments
		delta_nutrients = ((1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS) * self.fba.getExternalExchangeFluxes()).asNumber().astype(int)

		# TODO (Eran) use environment_molecule_ids rather than external_exchange_molecules, delta_nutrients needs to be for all env molecules
		self.environment_molecules.countsInc(self.external_exchange_molecule_ids, delta_nutrients)

		# Write outputs to listeners
		self.writeToListener("FBAResults", "deltaMetabolites", metaboliteCountsFinal - metaboliteCountsInit)
		self.writeToListener("FBAResults", "reactionFluxes", self.fba.getReactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes", exFluxes)
		self.writeToListener("FBAResults", "objectiveValue", self.fba.getObjectiveValue())
		self.writeToListener("FBAResults", "shadowPrices", self.fba.getShadowPrices(self.metaboliteNames))
		self.writeToListener("FBAResults", "reducedCosts", self.fba.getReducedCosts(self.fba.getReactionIDs()))
		self.writeToListener("FBAResults", "targetConcentrations", [self.homeostaticObjective[mol] for mol in self.fba.getHomeostaticTargetMolecules()])
		self.writeToListener("FBAResults", "homeostaticObjectiveValues", self.fba.getHomeostaticObjectiveValues())
		self.writeToListener("FBAResults", "kineticObjectiveValues", self.fba.getKineticObjectiveValues())

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit", metaboliteCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal", metaboliteCountsFinal)
		self.writeToListener("EnzymeKinetics", "enzymeCountsInit", kineticsEnzymesCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations", metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		self.writeToListener("EnzymeKinetics", "countsToMolar", countsToMolar.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		self.writeToListener("EnzymeKinetics", "actualFluxes", self.fba.getReactionFluxes(self.kineticsConstrainedReactions) / self.timeStepSec())

		self.writeToListener("EnzymeKinetics", "targetFluxes", targets / self.timeStepSec())
		self.writeToListener("EnzymeKinetics", "reactionConstraint", reactionConstraint[self.active_constraints_mask])

	# limit amino acid uptake to what is needed to meet concentration objective to prevent use as carbon source
	def _setExternalMoleculeLevels(self, externalMoleculeLevels, metaboliteConcentrations):
		for aa in self.AAs:
			if aa + "[p]" in self.fba.getExternalMoleculeIDs():
				idx = self.externalMoleculeIDs.index(aa + "[p]")
			elif aa + "[c]" in self.fba.getExternalMoleculeIDs():
				idx = self.externalMoleculeIDs.index(aa + "[c]")
			else:
				continue

			concDiff = self.homeostaticObjective[aa + "[c]"] - metaboliteConcentrations[self.metaboliteNames.index(aa + "[c]")].asNumber(COUNTS_UNITS / VOLUME_UNITS)
			if concDiff < 0:
				concDiff = 0

			if externalMoleculeLevels[idx] > concDiff:
				externalMoleculeLevels[idx] =  concDiff

		self.fba.setExternalMoleculeLevels(externalMoleculeLevels)
