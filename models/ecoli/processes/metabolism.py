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

		# Local sim_data references
		metabolism = sim_data.process.metabolism
		constants = sim_data.constants
		mass = sim_data.mass

		# Load constants
		self.nAvogadro = constants.nAvogadro
		self.cellDensity = constants.cellDensity
		self.ngam = constants.nonGrowthAssociatedMaintenance
		energyCostPerWetMass = constants.darkATP * mass.cellDryMassFraction

		self.exchangeConstraints = metabolism.exchangeConstraints

		self._biomass_concentrations = {}
		self._getBiomassAsConcentrations = mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		self.use_trna_charging = sim._trna_charging

		# Include ppGpp concentration target in objective if not handled kinetically in other processes
		self.include_ppgpp = not sim._ppgpp_regulation or not self.use_trna_charging
		self.ppgpp_id = sim_data.moleculeIds.ppGpp
		self.getppGppConc = sim_data.growthRateParameters.getppGppConc

		# Use information from the environment
		environment = self._external_states['Environment']

		# go through all media in the timeline and add to metaboliteNames
		self.metaboliteNamesFromNutrients = set()
		exchange_molecules = set()
		if self.include_ppgpp:
			self.metaboliteNamesFromNutrients.add(self.ppgpp_id)
		for time, media_id in environment.current_timeline:
			self.metaboliteNamesFromNutrients.update(
				metabolism.concentrationUpdates.concentrationsBasedOnNutrients(media_id)
				)
			exchanges = sim_data.external_state.exchange_data_from_media(media_id)
			exchange_molecules.update(exchanges['externalExchangeMolecules'])
		self.metaboliteNamesFromNutrients = list(sorted(self.metaboliteNamesFromNutrients))
		exchange_molecules = list(sorted(exchange_molecules))
		moleculeMasses = dict(zip(exchange_molecules,
			sim_data.getter.getMass(exchange_molecules).asNumber(MASS_UNITS / COUNTS_UNITS)))

		concDict = metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			environment.current_media_id
			)
		doubling_time = sim_data.conditionToDoublingTime[sim_data.condition]
		concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			doubling_time
			)
		concDict.update(concModificationsBasedOnCondition)

		if self.include_ppgpp:
			concDict[self.ppgpp_id] = self.getppGppConc(doubling_time)

		self.homeostaticObjective = dict((key, concDict[key].asNumber(CONC_UNITS)) for key in concDict)


		# Data structures to compute reaction bounds based on enzyme presence/absence
		self.catalyst_ids = metabolism.catalyst_ids
		self.reactions_with_catalyst = metabolism.reactions_with_catalyst

		catalysisMatrixI = metabolism.catalysisMatrixI
		catalysisMatrixJ = metabolism.catalysisMatrixJ
		catalysisMatrixV = metabolism.catalysisMatrixV

		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
		self.catalysisMatrix = csr_matrix((catalysisMatrixV, (catalysisMatrixI, catalysisMatrixJ)), shape = shape)

		# Function to compute reaction targets based on kinetic parameters and molecule concentrations
		self.getKineticConstraints = metabolism.getKineticConstraints

		# Remove disabled reactions so they don't get included in the FBA problem setup
		if hasattr(metabolism, "kineticTargetShuffleRxns") and metabolism.kineticTargetShuffleRxns is not None:
			self.kinetics_constrained_reactions = metabolism.kineticTargetShuffleRxns
			self.active_constraints_mask = np.ones(len(self.kinetics_constrained_reactions), dtype=bool)
		else:
			kinetic_constraint_reactions = metabolism.kinetic_constraint_reactions
			constraintsToDisable = metabolism.constraintsToDisable
			self.active_constraints_mask = np.array([(rxn not in constraintsToDisable) for rxn in kinetic_constraint_reactions])
			self.kinetics_constrained_reactions = list(np.array(kinetic_constraint_reactions)[self.active_constraints_mask])

		# Add kinetic reaction targets from boundary
		self.boundary_constrained_reactions = environment.transport_fluxes.keys()
		self.all_constrained_reactions = self.kinetics_constrained_reactions + self.boundary_constrained_reactions

		self.kinetic_constraint_enzymes = metabolism.kinetic_constraint_enzymes
		self.kinetic_constraint_substrates = metabolism.kinetic_constraint_substrates

		# Set solver and kinetic objective weight (lambda)
		solver = metabolism.solver
		kinetic_objective_weight = metabolism.kinetic_objective_weight
		kinetic_objective_weight_in_range = metabolism.kinetic_objective_weight_in_range

		# Disable kinetics completely if weight is 0 or specified in file above
		if not USE_KINETICS or kinetic_objective_weight == 0:
			objective_type = 'homeostatic'
			self.use_kinetics = False
			kinetic_objective_weight = 0
		else:
			objective_type = 'homeostatic_kinetics_mixed'
			self.use_kinetics = True

		# Set up FBA solver
		# reactionRateTargets value is just for initialization, it gets reset each timestep during evolveState
		self.fbaObjectOptions = {
			"reactionStoich": metabolism.reactionStoich,
			"externalExchangedMolecules": exchange_molecules,
			"objective": self.homeostaticObjective,
			"objectiveType": objective_type,
			"objectiveParameters": {
					"kineticObjectiveWeight": kinetic_objective_weight,
					'kinetic_objective_weight_in_range': kinetic_objective_weight_in_range,
					"reactionRateTargets": {reaction: 1 for reaction in self.all_constrained_reactions},
					"oneSidedReactionTargets": [],
					},
			"moleculeMasses": moleculeMasses,
			"secretionPenaltyCoeff": metabolism.secretion_penalty_coeff, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver": solver,
			"maintenanceCostGAM": energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction": metabolism.maintenanceReaction,
		}
		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])

		# Disable all rates during burn-in
		if self.use_kinetics:
			if KINETICS_BURN_IN_PERIOD > 0:
				self.fba.disableKineticTargets()
				self.burnInComplete = False
			else:
				self.burnInComplete = True

		## Construct views
		# views on metabolism bulk molecules
		self.metaboliteNames = self.fba.getOutputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.catalyst_ids)
		self.kineticsEnzymes = self.bulkMoleculesView(self.kinetic_constraint_enzymes)
		self.kineticsSubstrates = self.bulkMoleculesView(self.kinetic_constraint_substrates)

		assert self.metaboliteNames == self.fba.getInternalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		self.aa_names_no_location = [x[:-3] for x in sorted(sim_data.amino_acid_1_to_3_ordered.values())]

		self.shuffleIdxs = None
		if hasattr(metabolism, "kineticTargetShuffleIdxs") and metabolism.kineticTargetShuffleIdxs is not None:
			self.shuffleIdxs = metabolism.kineticTargetShuffleIdxs

		self.shuffleCatalyzedIdxs = None
		if hasattr(metabolism, "catalystShuffleIdxs") and metabolism.catalystShuffleIdxs is not None:
			self.shuffleCatalyzedIdxs = metabolism.catalystShuffleIdxs

		# Track updated AA concentration targets with tRNA charging
		self.aa_targets = {}
		self.aa_targets_not_updated = set(['L-SELENOCYSTEINE[c]'])
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

		# Get environment updates
		environment = self._external_states['Environment']
		current_media_id = environment.current_media_id
		exchange_data = environment.get_exchange_data()

		# make sure there are no new flux targets from the boundary
		assert set(environment.transport_fluxes.keys()).issubset(self.all_constrained_reactions)

		doubling_time = self.nutrientToDoublingTime.get(current_media_id, self.nutrientToDoublingTime["minimal"])
		concModificationsBasedOnCondition = self.getBiomassAsConcentrations(doubling_time)

		if self.use_trna_charging:
			concModificationsBasedOnCondition.update(self.update_amino_acid_targets(countsToMolar))
		if self.include_ppgpp:
			concModificationsBasedOnCondition[self.ppgpp_id] = self.getppGppConc(doubling_time)

		# Set external molecule levels
		external_exchange_molecule_ids = self.fba.getExternalMoleculeIDs()
		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			external_exchange_molecule_ids,
			coefficient,
			CONC_UNITS,
			current_media_id,
			exchange_data,
			concModificationsBasedOnCondition,
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

		# Maintenance reactions
		## Calculate new NGAM
		flux = (self.ngam * coefficient).asNumber(CONC_UNITS)
		self.fba.setReactionFluxBounds(
			self.fba._reactionID_NGAM,
			lowerBounds=flux, upperBounds=flux,
			)

		## Calculate GTP usage based on how much was needed in polypeptide
		## elongation in previous step.
		gtp_to_hydrolyze = self._sim.processes["PolypeptideElongation"].gtp_to_hydrolyze
		flux = (countsToMolar * gtp_to_hydrolyze).asNumber(CONC_UNITS)
		self.fba.setReactionFluxBounds(
			self.fba._reactionID_polypeptideElongationEnergy,
			lowerBounds=flux, upperBounds=flux,
			)

		# Constrain reactions based on absence of catalysts
		## Read counts for catalysts and enzymes (catalysts with kinetics constraints)
		catalystsCountsInit = self.catalysts.counts()

		## Set hard upper bounds constraints based on enzyme presence (infinite upper bound) or absence (upper bound of zero)
		catalyzedReactionBounds = np.inf * np.ones(len(self.reactions_with_catalyst))
		rxnPresence = self.catalysisMatrix.dot(catalystsCountsInit)
		catalyzedReactionBounds[rxnPresence == 0] = 0
		if self.shuffleCatalyzedIdxs is not None:
			catalyzedReactionBounds = catalyzedReactionBounds[self.shuffleCatalyzedIdxs]
		self.fba.setReactionFluxBounds(self.reactions_with_catalyst,
			upperBounds=catalyzedReactionBounds, raiseForReversible=False)

		# Constrain reactions based on kinetic values
		kineticsEnzymesCountsInit = self.kineticsEnzymes.counts()
		kineticsEnzymesConcentrations = countsToMolar * kineticsEnzymesCountsInit

		kineticsSubstratesCountsInit = self.kineticsSubstrates.counts()
		kineticsSubstratesConcentrations = countsToMolar * kineticsSubstratesCountsInit

		## Set target fluxes for reactions based on their most relaxed constraint
		reactionTargets = (units.umol / units.L / units.s) * self.getKineticConstraints(
			kineticsEnzymesConcentrations.asNumber(units.umol / units.L),
			kineticsSubstratesConcentrations.asNumber(units.umol / units.L),
			)

		## Shuffle parameters (only performed in very specific cases)
		if self.shuffleIdxs is not None:
			reactionTargets = (units.umol / units.L / units.s) * reactionTargets.asNumber()[self.shuffleIdxs, :]

		## Calculate reaction flux target for current time step
		targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(CONC_UNITS)[self.active_constraints_mask, :]

		# add boundary targets
		transport_targets = environment.transport_fluxes.values()
		lower_targets = np.concatenate((targets[:, 0], transport_targets), axis=0)
		mean_targets = np.concatenate((targets[:, 1], transport_targets), axis=0)
		upper_targets = np.concatenate((targets[:, 2], transport_targets), axis=0)

		## Set kinetic targets only if kinetics is enabled
		if self.use_kinetics and self.burnInComplete:
			self.fba.setKineticTarget(
				self.all_constrained_reactions, mean_targets,
				lower_targets=lower_targets, upper_targets=upper_targets)

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
		environment.molecule_exchange(external_exchange_molecule_ids, delta_nutrients)

		import_exchange, import_constraint = environment.get_import_constraints(exchange_data)

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
		self.writeToListener("EnzymeKinetics", "targetFluxes", mean_targets / self.timeStepSec())
		# TODO: add lower and upper targets

	# limit amino acid uptake to what is needed to meet concentration objective to prevent use as carbon source
	def _setExternalMoleculeLevels(self, externalMoleculeLevels, metaboliteConcentrations):
		external_exchange_molecule_ids = self.fba.getExternalMoleculeIDs()
		for aa in self.aa_names_no_location:
			if aa + "[p]" in external_exchange_molecule_ids:
				idx = external_exchange_molecule_ids.index(aa + "[p]")
			elif aa + "[c]" in external_exchange_molecule_ids:
				idx = external_exchange_molecule_ids.index(aa + "[c]")
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

		if doubling_time not in self._biomass_concentrations:
			self._biomass_concentrations[doubling_time] = self._getBiomassAsConcentrations(doubling_time)

		return self._biomass_concentrations[doubling_time]

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

		TODO:
		- remove access to PolypeptideElongation class attribute (aa_count_diff)
		'''

		count_diff = self._sim.processes['PolypeptideElongation'].aa_count_diff

		if len(self.aa_targets):
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
