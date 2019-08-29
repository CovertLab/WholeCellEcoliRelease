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
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

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

		self.biomass_concentrations = {}
		self._getBiomassAsConcentrations = sim_data.mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		# create boundary object
		self.boundary = Boundary(
			sim_data.external_state.environment,
			sim_data.process.metabolism.boundary,
			self._external_states,
			self.environmentView
			)

		# go through all media in the timeline and add to metaboliteNames
		self.metaboliteNamesFromNutrients = set()
		for time, media_id in self.boundary.current_timeline:
			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					media_id, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			self.boundary.current_media
			)
		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			sim_data.conditionToDoublingTime[sim_data.condition]
			)
		concDict.update(self.concModificationsBasedOnCondition)
		self.homeostaticObjective = dict((key, concDict[key].asNumber(CONC_UNITS)) for key in concDict)

		# Load initial mass
		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit
		initCellMass = initWaterMass + initDryMass
		energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass
		moleculeMasses = dict(zip(self.boundary.exchange_data['externalExchangeMolecules'],
			sim_data.getter.getMass(self.boundary.exchange_data['externalExchangeMolecules']).asNumber(MASS_UNITS / COUNTS_UNITS)))

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
			self.kinetics_constrained_reactions = sim_data.process.metabolism.kineticTargetShuffleRxns
			self.active_constraints_mask = np.ones(len(self.kinetics_constrained_reactions), dtype=bool)
		else:
			constrainedReactionList = sim_data.process.metabolism.constrainedReactionList
			constraintsToDisable = sim_data.process.metabolism.constraintsToDisable
			self.active_constraints_mask = np.array([(rxn not in constraintsToDisable) for rxn in constrainedReactionList])
			self.kinetics_constrained_reactions = list(np.array(constrainedReactionList)[self.active_constraints_mask])

		# Add kinetic reaction targets from boundary
		self.boundary_constrained_reactions = self.boundary.transport_fluxes.keys()
		self.all_constrained_reactions = self.kinetics_constrained_reactions + self.boundary_constrained_reactions

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
			"externalExchangedMolecules" : self.boundary.exchange_data['externalExchangeMolecules'],
			"objective" : self.homeostaticObjective,
			"objectiveType" : "homeostatic_kinetics_mixed",
			"objectiveParameters" : {
					"kineticObjectiveWeight" : kinetic_objective_weight,
					"reactionRateTargets" : {reaction: 1 for reaction in self.all_constrained_reactions},
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
		self.currentNgam = 1 * CONC_UNITS
		self.currentPolypeptideElongationEnergy = 1 * CONC_UNITS

		# External molecules
		self.externalMoleculeIDs = self.fba.getExternalMoleculeIDs()

		## Construct views
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

		# Track updated AA concentration targets with tRNA charging
		self.use_trna_charging = sim._trna_charging
		self.aa_targets = {}
		self.aa_targets_not_updated = set(['L-SELENOCYSTEINE[c]', 'GLT[c]'])
		self.aa_names = sim_data.moleculeGroups.aaIDs
		self.aas = self.bulkMoleculesView(self.aa_names)

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

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)

		# get boundary conditions
		self.boundary.updateBoundary()
		current_media = self.boundary.current_media
		exchange_data = self.boundary.exchange_data

		# make sure there are no new flux targets from the boundary
		assert set(self.boundary.transport_fluxes.keys()).issubset(self.all_constrained_reactions)

		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			self.nutrientToDoublingTime.get(current_media, self.nutrientToDoublingTime["minimal"])
			)

		if self.use_trna_charging:
			self.concModificationsBasedOnCondition.update(self.update_amino_acid_targets(countsToMolar))

		# Set external molecule levels
		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			CONC_UNITS,
			current_media,
			exchange_data,
			self.concModificationsBasedOnCondition,
			)

		if newObjective != None and newObjective != self.homeostaticObjective:
			self.fba.update_homeostatic_targets(newObjective)
			self.homeostaticObjective = newObjective

		# After completing the burn-in, enable kinetic rates
		if self.use_kinetics and (not self.burnInComplete) and (self._sim.time() > KINETICS_BURN_IN_PERIOD):
			self.burnInComplete = True
			self.fba.enableKineticTargets()

		#  Find metabolite concentrations from metabolite counts
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit[self.internalExchangeIdxs]

		# Make a dictionary of metabolite names to metabolite concentrations
		self.fba.setInternalMoleculeLevels(metaboliteConcentrations.asNumber(CONC_UNITS))

		# Set external molecule levels
		# TODO -- this can change external AA levels for the fba problem. Problematic for reliable control of environmental response
		self._setExternalMoleculeLevels(externalMoleculeLevels, metaboliteConcentrations)

		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		ADJUSTMENT_RATIO = .01

		# Calculate new NGAM and update if necessary
		self.newNgam = self.ngam * coefficient
		ngam_diff = np.abs((self.currentNgam - self.newNgam).asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		if ngam_diff > ADJUSTMENT_RATIO:
			self.currentNgam = self.newNgam
			flux = (self.ngam * coefficient).asNumber(CONC_UNITS)
			self.fba.setReactionFluxBounds(self.fba._reactionID_NGAM, lowerBounds=flux, upperBounds=flux)

		# Calculate GTP usage based on how much was needed in polypeptide elongation in previous step and update if necessary
		newPolypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			newPolypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest
		poly_diff = np.abs((self.currentPolypeptideElongationEnergy - newPolypeptideElongationEnergy).asNumber()) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		if poly_diff > ADJUSTMENT_RATIO:
			self.currentPolypeptideElongationEnergy = newPolypeptideElongationEnergy
			flux = self.currentPolypeptideElongationEnergy.asNumber(CONC_UNITS)
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
		targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(CONC_UNITS)[self.active_constraints_mask]

		# add boundary targets
		all_targets = np.concatenate((targets, self.boundary.transport_fluxes.values()), axis=0)

		## Set kinetic targets only if kinetics is enabled
		if self.use_kinetics and self.burnInComplete:
			self.fba.setKineticTarget(self.all_constrained_reactions, all_targets, raiseForReversible = False)

		# Solve FBA problem and update metabolite counts
		deltaMetabolites = (1 / countsToMolar) * (CONC_UNITS * self.fba.getOutputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		exchange_fluxes = CONC_UNITS * self.fba.getExternalExchangeFluxes()
		converted_exchange_fluxes = (exchange_fluxes / coefficient).asNumber(units.mmol / units.g / units.h)

		# update environmental nutrient counts
		delta_nutrients = ((1 / countsToMolar) * exchange_fluxes).asNumber().astype(int)
		external_exchange_molecule_ids = self.fba.getExternalMoleculeIDs()
		self.boundary.updateEnvironment(external_exchange_molecule_ids, delta_nutrients)

		import_exchange, import_constraint = self.boundary.getImportConstraints(exchange_data)

		# Write outputs to listeners
		self.writeToListener("FBAResults", "import_exchange", import_exchange)
		self.writeToListener("FBAResults", "import_constraint", import_constraint)
		self.writeToListener("FBAResults", "deltaMetabolites", metaboliteCountsFinal - metaboliteCountsInit)
		self.writeToListener("FBAResults", "reactionFluxes", self.fba.getReactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes", converted_exchange_fluxes)
		self.writeToListener("FBAResults", "objectiveValue", self.fba.getObjectiveValue())
		self.writeToListener("FBAResults", "shadowPrices", self.fba.getShadowPrices(self.metaboliteNames))
		self.writeToListener("FBAResults", "reducedCosts", self.fba.getReducedCosts(self.fba.getReactionIDs()))
		self.writeToListener("FBAResults", "targetConcentrations", [self.homeostaticObjective[mol] for mol in self.fba.getHomeostaticTargetMolecules()])
		self.writeToListener("FBAResults", "homeostaticObjectiveValues", self.fba.getHomeostaticObjectiveValues())
		self.writeToListener("FBAResults", "kineticObjectiveValues", self.fba.getKineticObjectiveValues())

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit", metaboliteCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal", metaboliteCountsFinal)
		self.writeToListener("EnzymeKinetics", "enzymeCountsInit", kineticsEnzymesCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations", metaboliteConcentrations.asNumber(CONC_UNITS))
		self.writeToListener("EnzymeKinetics", "countsToMolar", countsToMolar.asNumber(CONC_UNITS))
		self.writeToListener("EnzymeKinetics", "actualFluxes", self.fba.getReactionFluxes(self.all_constrained_reactions) / self.timeStepSec())
		self.writeToListener("EnzymeKinetics", "targetFluxes", all_targets / self.timeStepSec())
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

			concDiff = self.homeostaticObjective[aa + "[c]"] - metaboliteConcentrations[self.metaboliteNames.index(aa + "[c]")].asNumber(CONC_UNITS)
			if concDiff < 0:
				concDiff = 0

			if externalMoleculeLevels[idx] > concDiff:
				externalMoleculeLevels[idx] =  concDiff

		self.fba.setExternalMoleculeLevels(externalMoleculeLevels)

	def getBiomassAsConcentrations(self, doubling_time):
		'''
		Caches the result of the sim_data function to improve performance since
		function requires computation but won't change for a given doubling_time.

		Args:
			doubling_time (float with time units): doubling time of the cell to
				get the metabolite concentrations for

		Returns:
			dict {str : float with concentration units}: dictionary with metabolite
				IDs as keys and concentrations as values
		'''

		if doubling_time not in self.biomass_concentrations:
			self.biomass_concentrations[doubling_time] = self._getBiomassAsConcentrations(doubling_time)

		return self.biomass_concentrations[doubling_time]

	def update_amino_acid_targets(self, counts_to_molar):
		'''
		Finds new amino acid concentration targets based on difference in supply
		and number of amino acids used in polypeptide_elongation

		Args:
			counts_to_molar (float with mol/volume units): conversion from counts
				to molar for the current state of the cell

		Returns:
			dict {AA name (str): AA conc (float with mol/volume units)}:
				new concentration targets for each amino acid

		Skips updates to certain molecules defined in self.aa_targets_not_updated:
		- L-SELENOCYSTEINE: rare amino acid that led to high variability when updated
		- GLT: high measured concentration that never doubled causing slow growth
		- LEU: increase in concentration caused TF regulation to stop transcription
		  of AA synthesis pathway genes

		TODO:
		- remove access to PolypeptideElongation class attribute (aa_count_diff)
		'''

		count_diff = self._sim.processes['PolypeptideElongation'].aa_count_diff

		if len(count_diff):
			for aa, diff in count_diff.items():
				if aa in self.aa_targets_not_updated:
					continue
				self.aa_targets[aa] += diff
		# First time step of a simulation so set target to current counts to prevent
		# concentration jumps between generations
		else:
			for aa, counts in zip(self.aa_names, self.aas.total_counts()):
				if aa in self.aa_targets_not_updated:
					continue
				self.aa_targets[aa] = counts

		return {aa: counts * counts_to_molar for aa, counts in self.aa_targets.items()}


class Boundary(object):
	'''
	Boundary provides an interface between metabolism and the environment.
	This includes handling all references to the media condition and exchange_data
	'''

	def __init__(self, sim_data_environment, sim_data_boundary, external_state, environmentView):
		self.external_state = external_state

		# get maps between environment and exchange molecules
		self.env_to_exchange_map = sim_data_environment.env_to_exchange_map
		self.exchange_to_env_map = sim_data_environment.exchange_to_env_map

		# get functions
		self.exchangeDataFromConcentrations = sim_data_boundary.exchangeDataFromConcentrations
		self.exchangeDataFromMedia = sim_data_boundary.exchangeDataFromMedia
		self.getImportConstraints = sim_data_boundary.getImportConstraints

		# get variables from environment
		self.current_timeline = self.external_state['Environment'].current_timeline

		# views on environment
		self.environment_molecule_ids = self.external_state['Environment']._moleculeIDs
		self.environment_molecules = environmentView(self.environment_molecule_ids)

		# transport fluxes from the external state
		self.transport_fluxes = self.external_state['Environment'].transport_fluxes

		self.updateBoundary()

	def updateBoundary(self):
		'''
		update all boundary variables for the current environment
		'''

		self.current_media = self.external_state['Environment'].current_media_id
		current_concentrations = dict(zip(self.environment_molecule_ids, self.environment_molecules.totalConcentrations()))
		self.exchange_data = self.exchangeDataFromConcentrations(current_concentrations)

		# transport fluxes from the external state
		self.transport_fluxes = self.external_state['Environment'].transport_fluxes

	def updateEnvironment(self, external_exchange_molecule_ids, delta_nutrients):
		'''
		Convert exchange molecules to environmental molecules using mapping, and passes delta counts to the local environment

		Args:
			external_exchange_molecule_ids (tuple[str]): a tuple with all the external exchange molecules from FBA
			delta_nutrients (np.array[int]): an array with the delta counts for all of the exchange molecules.
				The array length should be the same as external_exchange_molecule_ids.

		Modifies:
			uses countsInc() to update local_environment's _env_delta_counts
		'''

		mapped_environment_molecule_ids = [self.exchange_to_env_map[mol_id] for mol_id in external_exchange_molecule_ids]
		self.environment_molecules.countsInc(mapped_environment_molecule_ids, delta_nutrients)
