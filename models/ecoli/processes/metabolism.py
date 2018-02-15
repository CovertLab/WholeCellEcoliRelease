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
from wholecell.utils.enzymeKinetics import EnzymeKinetics

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

SECRETION_PENALTY_COEFF = 1e-5

NONZERO_ENZYMES = False

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

		self.nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel

		self.getBiomassAsConcentrations = sim_data.mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		# Create objective for homeostatic constraints
		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel][0][1]
			)
		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[sim_data.condition])
		concDict.update(self.concModificationsBasedOnCondition)
		self.homeostaticObjective = dict((key, concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in concDict)

		# Load initial mass
		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit
		initCellMass = initWaterMass + initDryMass

		energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		# Setup molecules in external environment that can be exchanged
		externalExchangedMolecules = sim_data.nutrientData["secretionExchangeMolecules"]
		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel]:
			externalExchangedMolecules += sim_data.nutrientData["importExchangeMolecules"][nutrientsLabel]

			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					nutrientsLabel, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		externalExchangedMolecules = sorted(externalExchangedMolecules)
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

		# Data structures to compute reaction targets based on kinetic parameters
		from reconstruction.ecoli.dataclasses.process.metabolism_constraints import constraints as kineticsConstraints
		self.kineticsConstraints = kineticsConstraints

		self.useAllConstraints = sim_data.process.metabolism.useAllConstraints
		self.constraintsToDisable = sim_data.process.metabolism.constraintsToDisable
		self.kineticsConstrainedReactions = sim_data.process.metabolism.constrainedReactionList
		if hasattr(sim_data.process.metabolism, "kineticTargetShuffleRxns") and sim_data.process.metabolism.kineticTargetShuffleRxns != None:
			self.kineticsConstrainedReactions = sim_data.process.metabolism.kineticTargetShuffleRxns
			self.useAllConstraints = True
		self.kineticsEnzymesList = sim_data.process.metabolism.enzymeIdList
		self.kineticsSubstratesList = sim_data.process.metabolism.kineticsSubstratesList

		constraintToReactionMatrixI = sim_data.process.metabolism.constraintToReactionMatrixI
		constraintToReactionMatrixJ = sim_data.process.metabolism.constraintToReactionMatrixJ
		constraintToReactionMatrixV = sim_data.process.metabolism.constraintToReactionMatrixV
		shape = (constraintToReactionMatrixI.max() + 1, constraintToReactionMatrixJ.max() + 1)
		self.constraintToReactionMatrix = np.zeros(shape, np.float64)
		self.constraintToReactionMatrix[constraintToReactionMatrixI, constraintToReactionMatrixJ] = constraintToReactionMatrixV
		self.constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly

		# Set up FBA solver
		# reactionRateTargets value is just for initialization, it gets reset each timestep during evolveState
		self.fbaObjectOptions = {
			"reactionStoich" : sim_data.process.metabolism.reactionStoich,
			"externalExchangedMolecules" : externalExchangedMolecules,
			"objective" : self.homeostaticObjective,
			"objectiveType" : "homeostatic_kinetics_mixed",
			"objectiveParameters" : {
					"kineticObjectiveWeight" : sim_data.constants.metabolismKineticObjectiveWeight,
					"reactionRateTargets" : {reaction : 1 for reaction in self.kineticsConstrainedReactions},
					"oneSidedReactionTargets" : [],
					},
			"moleculeMasses" : moleculeMasses,
			"secretionPenaltyCoeff" : SECRETION_PENALTY_COEFF, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver" : "glpk",
			"maintenanceCostGAM" : energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction" : sim_data.process.metabolism.maintenanceReaction,
		}
		if USE_KINETICS == False:
			self.fbaObjectOptions["objectiveType"] = "homeostatic"
		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.outputMoleculeIDs()])

		# Disable all rates during burn-in
		if USE_KINETICS:
			if KINETICS_BURN_IN_PERIOD > 0:
				self.fba.disableKineticTargets()
				self.burnInComplete = False
			else:
				self.burnInComplete = True
				if not self.useAllConstraints:
					for rxn in self.constraintsToDisable:
						self.fba.disableKineticTargets(rxn)

		# Values will get updated at each time point
		self.currentNgam = 1 * (COUNTS_UNITS / VOLUME_UNITS)
		self.currentPolypeptideElongationEnergy = 1 * (COUNTS_UNITS / VOLUME_UNITS)

		# External molecules
		self.externalMoleculeIDs = self.fba.externalMoleculeIDs()

		# Views
		self.metaboliteNames = self.fba.outputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.catalystsList)
		self.kineticsEnzymes = self.bulkMoleculesView(self.kineticsEnzymesList)
		self.kineticsSubstrates = self.bulkMoleculesView(self.kineticsSubstratesList)

		outputMoleculeIDs = self.fba.outputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.internalMoleculeIDs()

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

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(self.nutrientToDoublingTime.get(self._sim.processes["PolypeptideElongation"].currentNutrients, self.nutrientToDoublingTime["minimal"]))

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)

		# Set external molecule levels
		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			self.nutrientsTimeSeriesLabel,
			self.time(),
			self.concModificationsBasedOnCondition,
			)

		updatedObjective = False
		if newObjective != None and newObjective != self.homeostaticObjective:
			# Build new fba instance with new objective
			self.fbaObjectOptions["objective"] = newObjective
			self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
			self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.outputMoleculeIDs()])
			self.homeostaticObjective = newObjective
			updatedObjective = True

		# After completing the burn-in, enable kinetic rates
		if (USE_KINETICS) and (not self.burnInComplete) and (self._sim.time() > KINETICS_BURN_IN_PERIOD):
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
		self.fba.internalMoleculeLevelsIs(metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# Set external molecule levels
		self._setExternalMoleculeLevels(externalMoleculeLevels, metaboliteConcentrations)

		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		ADJUSTMENT_RATIO = .01

		# Calculate new NGAM and update if necessary
		self.newNgam = self.ngam * coefficient
		ngam_diff = np.abs(self.currentNgam.asNumber() - self.newNgam.asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		if ngam_diff > ADJUSTMENT_RATIO:
			self.currentNgam = self.newNgam
			self.fba.maxReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.fba.minReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# Calculate GTP usage based on how much was needed in polypeptide elongation in previous step and update if necessary
		newPolypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			newPolypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest
		poly_diff = np.abs((self.currentPolypeptideElongationEnergy.asNumber() - newPolypeptideElongationEnergy.asNumber())) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		if poly_diff > ADJUSTMENT_RATIO:
			self.currentPolypeptideElongationEnergy = newPolypeptideElongationEnergy
			self.fba.maxReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))
			self.fba.minReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# Read counts for catalysts and enzymes (catalysts with kinetics constraints)
		catalystsCountsInit = self.catalysts.counts()

		kineticsEnzymesCountsInit = self.kineticsEnzymes.counts()
		kineticsEnzymesConcentrations = countsToMolar * kineticsEnzymesCountsInit

		kineticsSubstratesCountsInit = self.kineticsSubstrates.counts()
		kineticsSubstratesConcentrations = countsToMolar * kineticsSubstratesCountsInit

		# Add one of every enzyme to ensure none are zero
		if NONZERO_ENZYMES:
			catalystsCountsInit += 1
			kineticsEnzymesConcentrations = countsToMolar * (kineticsEnzymesCountsInit + 1)

		# Calculate and set kinetic targets if kinetics is enabled
		if USE_KINETICS and self.burnInComplete:
			# Set hard upper bounds constraints based on enzyme presence (infinite upper bound) or absence (upper bound of zero)
			catalyzedReactionBounds = np.inf * np.ones(len(self.reactionsWithCatalystsList))
			rxnPresence = self.catalysisMatrix.dot(catalystsCountsInit)
			catalyzedReactionBounds[rxnPresence == 0] = 0
			if self.shuffleCatalyzedIdxs is not None:
				catalyzedReactionBounds = catalyzedReactionBounds[self.shuffleCatalyzedIdxs]


			# Only update reaction limits that are different from previous time step
			updateIdxs = np.where(catalyzedReactionBounds != self.catalyzedReactionBoundsPrev)[0]
			updateRxns = [self.reactionsWithCatalystsList[idx] for idx in updateIdxs]
			updateVals = catalyzedReactionBounds[updateIdxs]
			self.fba.setMaxReactionFluxes(updateRxns, updateVals, raiseForReversible = False)
			self.catalyzedReactionBoundsPrev = catalyzedReactionBounds

			# Set target fluxes for reactions based on their most relaxed constraint
			# kineticsConstraints function created based on units of umol/L/s
			constraintValues = self.kineticsConstraints(
				kineticsEnzymesConcentrations.asNumber(units.umol / units.L),
				kineticsSubstratesConcentrations.asNumber(units.umol / units.L),
				)
			reactionTargets = (units.umol / units.L / units.s) * np.max(self.constraintToReactionMatrix * constraintValues, axis = 1)

			# Shuffle parameters (only performed in very specific cases)
			if self.shuffleIdxs is not None:
				reactionTargets = (units.umol / units.L / units.s) * reactionTargets.asNumber()[self.shuffleIdxs]

			# record which constraint was used, add constraintToReactionMatrix to ensure the index is one of the constraints if multiplication is 0
			reactionConstraint = np.argmax(self.constraintToReactionMatrix * constraintValues + self.constraintToReactionMatrix, axis = 1)

			# set reaction flux targets
			targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(COUNTS_UNITS / VOLUME_UNITS)
			self.fba.setKineticTarget(self.kineticsConstrainedReactions, targets, raiseForReversible = False)

			if not self.useAllConstraints:
				self.fba.disableKineticTargets(self.constraintsToDisable)

		# Solve FBA problem and update metabolite counts
		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.outputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)
		if USE_KINETICS and self.burnInComplete:
			relError = np.abs((self.fba.reactionFluxes(self.kineticsConstrainedReactions) - targets) / (targets + 1e-15))

		exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.externalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)

		# Write outputs to listeners
		self.writeToListener("FBAResults", "deltaMetabolites", metaboliteCountsFinal - metaboliteCountsInit)
		self.writeToListener("FBAResults", "reactionFluxes", self.fba.reactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes", exFluxes)
		self.writeToListener("FBAResults", "objectiveValue", self.fba.objectiveValue())
		self.writeToListener("FBAResults", "rowDualValues", self.fba.rowDualValues(self.metaboliteNames))
		self.writeToListener("FBAResults", "columnDualValues", self.fba.columnDualValues(self.fba.reactionIDs()))
		self.writeToListener("FBAResults", "targetConcentrations", [self.homeostaticObjective[mol] for mol in self.fba.homeostaticTargetMolecules()])

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit", metaboliteCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal", metaboliteCountsFinal)
		self.writeToListener("EnzymeKinetics", "enzymeCountsInit", kineticsEnzymesCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations", metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		self.writeToListener("EnzymeKinetics", "countsToMolar", countsToMolar.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		self.writeToListener("EnzymeKinetics", "actualFluxes", self.fba.reactionFluxes(self.kineticsConstrainedReactions) / self.timeStepSec())

		if USE_KINETICS and self.burnInComplete:
			self.writeToListener("EnzymeKinetics", "targetFluxes", targets / self.timeStepSec())
			self.writeToListener("EnzymeKinetics", "reactionConstraint", reactionConstraint)

	# limit amino acid uptake to what is needed to meet concentration objective to prevent use as carbon source
	def _setExternalMoleculeLevels(self, externalMoleculeLevels, metaboliteConcentrations):
		for aa in self.AAs:
			if aa + "[p]" in self.fba.externalMoleculeIDs():
				idx = self.externalMoleculeIDs.index(aa + "[p]")
			elif aa + "[c]" in self.fba.externalMoleculeIDs():
				idx = self.externalMoleculeIDs.index(aa + "[c]")
			else:
				continue

			concDiff = self.homeostaticObjective[aa + "[c]"] - metaboliteConcentrations[self.metaboliteNames.index(aa + "[c]")].asNumber(COUNTS_UNITS / VOLUME_UNITS)
			if concDiff < 0:
				concDiff = 0

			if externalMoleculeLevels[idx] > concDiff:
				externalMoleculeLevels[idx] =  concDiff

		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)
