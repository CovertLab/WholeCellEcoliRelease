"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- option to call a reduced form of metabolism (assume optimal)
- handle oneSidedReaction constraints
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.sparse import csr_matrix
from typing import List, Optional, Set, Tuple

from reconstruction.ecoli.simulation_data import SimulationDataEcoli
import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM
from wholecell.utils.modular_fba import FluxBalanceAnalysis
from six.moves import zip


COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS
CONVERSION_UNITS = MASS_UNITS * TIME_UNITS / VOLUME_UNITS
GDCW_BASIS = units.mmol / units.g / units.h

USE_KINETICS = True


class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Metabolism, self).initialize(sim, sim_data)

		# Use information from the environment and sim
		self.get_import_constraints = sim_data.external_state.get_import_constraints
		self.nutrientToDoublingTime = sim_data.nutrient_to_doubling_time
		environment = self._external_states['Environment']
		self.use_trna_charging = sim._trna_charging
		self.include_ppgpp = not sim._ppgpp_regulation or not self.use_trna_charging
		self.mechanistic_aa_uptake = sim._mechanistic_aa_uptake
		metabolism = sim_data.process.metabolism

		# Create model to use to solve metabolism updates
		self.model = FluxBalanceAnalysisModel(
			sim_data,
			imports=environment.get_import_molecules(),
			timeline=environment.current_timeline,
			include_ppgpp=self.include_ppgpp,
			)

		# Save constants
		constants = sim_data.constants
		self.nAvogadro = constants.n_avogadro
		self.cellDensity = constants.cell_density

		# Track updated AA concentration targets with tRNA charging
		self.aa_targets = {}
		self.aa_targets_not_updated = {'L-SELENOCYSTEINE[c]'}
		self.aa_names = sim_data.molecule_groups.amino_acids

		# Construct views
		self.metabolites = self.bulkMoleculesView(self.model.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.model.catalyst_ids)
		self.kineticsEnzymes = self.bulkMoleculesView(self.model.kinetic_constraint_enzymes)
		self.kineticsSubstrates = self.bulkMoleculesView(self.model.kinetic_constraint_substrates)
		self.aas = self.bulkMoleculesView(self.aa_names)

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		# Molecules with concentration updates for listener
		self.linked_metabolites = metabolism.concentration_updates.linked_metabolites
		doubling_time = self.nutrientToDoublingTime.get(
			environment.current_media_id,
			self.nutrientToDoublingTime["minimal"])
		update_molecules = list(self.model.getBiomassAsConcentrations(doubling_time).keys())
		if self.use_trna_charging:
			update_molecules += [aa for aa in self.aa_names if aa not in self.aa_targets_not_updated]
			update_molecules += list(self.linked_metabolites.keys())
		if self.include_ppgpp:
			update_molecules += [self.model.ppgpp_id]
		self.conc_update_molecules = sorted(update_molecules)

		# Amino acids in media for import rates
		self.aa_exchange_names = np.array([
			sim_data.external_state.env_to_exchange_map[aa[:-3]]
			for aa in self.aa_names
			])
		self.aa_environment = self.environmentView([aa[:-3] for aa in self.aa_exchange_names])

		self.removed_aa_uptake = np.array([
			aa in self.aa_targets_not_updated
			for aa in self.aa_exchange_names
			])

		self.amino_acid_import = metabolism.amino_acid_import
		self.aa_transporters_names = metabolism.aa_transporters_names
		self.aa_transporters_container = self.bulkMoleculesView(self.aa_transporters_names)

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.catalysts.requestAll()
		self.aa_transporters_container.requestAll()
		self.kineticsEnzymes.requestAll()
		self.kineticsSubstrates.requestAll()

	def evolveState(self):
		# Load current state of the sim
		## Get internal state variables
		time_step = self.timeStepSec() * units.s
		metabolite_counts_init = self.metabolites.counts()
		catalyst_counts = self.catalysts.counts()
		kinetic_enzyme_counts = self.kineticsEnzymes.counts()
		kinetic_substrate_counts = self.kineticsSubstrates.counts()
		translation_gtp = self._sim.processes["PolypeptideElongation"].gtp_to_hydrolyze
		cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
		dry_mass = self.readFromListener("Mass", "dryMass") * units.fg

		## Get environment updates
		environment = self._external_states['Environment']
		current_media_id = environment.current_media_id
		unconstrained, constrained = environment.get_exchange_data()

		## Calculate state values
		cellVolume = cell_mass / self.cellDensity
		counts_to_molar = (1 / (self.nAvogadro * cellVolume)).asUnit(CONC_UNITS)

		## Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		coefficient = dry_mass / cell_mass * self.cellDensity * time_step

		## Determine updates to concentrations depending on the current state
		doubling_time = self.nutrientToDoublingTime.get(current_media_id, self.nutrientToDoublingTime["minimal"])
		conc_updates = self.model.getBiomassAsConcentrations(doubling_time)
		if self.use_trna_charging:
			conc_updates.update(self.update_amino_acid_targets(counts_to_molar))
		if self.include_ppgpp:
			conc_updates[self.model.ppgpp_id] = self.model.getppGppConc(doubling_time).asUnit(CONC_UNITS)
		## Converted from units to make reproduction from listener data
		## accurate to model results (otherwise can have floating point diffs)
		conc_updates = {
			met: conc.asNumber(CONC_UNITS)
			for met, conc in conc_updates.items()
			}

		aa_uptake_package = None
		if self.mechanistic_aa_uptake:
			aa_in_media = self.aa_environment.import_present()
			aa_in_media[self.removed_aa_uptake] = False
			import_rates = (counts_to_molar * self.timeStepSec() * self.amino_acid_import(
				aa_in_media, dry_mass, self.aa_transporters_container.counts(),
				self.mechanistic_aa_uptake)).asNumber(CONC_UNITS)
			aa_uptake_package=(import_rates[aa_in_media], self.aa_exchange_names[aa_in_media], True)

		# Update FBA problem based on current state
		## Set molecule availability (internal and external)
		self.model.set_molecule_levels(metabolite_counts_init, counts_to_molar,
			coefficient, current_media_id, unconstrained, constrained, conc_updates, aa_uptake_package)

		## Set reaction limits for maintenance and catalysts present
		self.model.set_reaction_bounds(catalyst_counts, counts_to_molar,
			coefficient, translation_gtp)

		## Constrain reactions based on targets
		targets, upper_targets, lower_targets = self.model.set_reaction_targets(kinetic_enzyme_counts,
			kinetic_substrate_counts, counts_to_molar, time_step)

		# Solve FBA problem and update states
		n_retries = 3
		fba = self.model.fba
		fba.solve(n_retries)

		## Internal molecule changes
		delta_metabolites = (1 / counts_to_molar) * (CONC_UNITS * fba.getOutputMoleculeLevelsChange())
		metabolite_counts_final = np.fmax(stochasticRound(
			self.randomState,
			metabolite_counts_init + delta_metabolites.asNumber()
			), 0).astype(np.int64)
		self.metabolites.countsIs(metabolite_counts_final)

		## Environmental changes
		exchange_fluxes = CONC_UNITS * fba.getExternalExchangeFluxes()
		converted_exchange_fluxes = (exchange_fluxes / coefficient).asNumber(GDCW_BASIS)
		delta_nutrients = ((1 / counts_to_molar) * exchange_fluxes).asNumber().astype(int)
		environment.molecule_exchange(fba.getExternalMoleculeIDs(), delta_nutrients)

		# Write outputs to listeners
		unconstrained, constrained, uptake_constraints = self.get_import_constraints(
			unconstrained, constrained, GDCW_BASIS)
		time_step_unitless = time_step.asNumber(TIME_UNITS)
		self.writeToListener('FBAResults', 'media_id', current_media_id)
		self.writeToListener('FBAResults', 'conc_updates',
			[conc_updates[m] for m in self.conc_update_molecules])
		self.writeToListener('FBAResults', 'catalyst_counts', catalyst_counts)
		self.writeToListener('FBAResults', 'translation_gtp', translation_gtp)
		self.writeToListener('FBAResults', 'coefficient', coefficient.asNumber(CONVERSION_UNITS))
		self.writeToListener('FBAResults', 'unconstrained_molecules', unconstrained)
		self.writeToListener('FBAResults', 'constrained_molecules', constrained)
		self.writeToListener('FBAResults', 'uptake_constraints', uptake_constraints)
		self.writeToListener("FBAResults", "deltaMetabolites", metabolite_counts_final - metabolite_counts_init)
		self.writeToListener("FBAResults", "reactionFluxes", fba.getReactionFluxes() / time_step_unitless)
		self.writeToListener("FBAResults", "externalExchangeFluxes", converted_exchange_fluxes)
		self.writeToListener("FBAResults", "objectiveValue", fba.getObjectiveValue())
		self.writeToListener("FBAResults", "shadowPrices", fba.getShadowPrices(self.model.metaboliteNamesFromNutrients))
		self.writeToListener("FBAResults", "reducedCosts", fba.getReducedCosts(fba.getReactionIDs()))
		self.writeToListener("FBAResults", "targetConcentrations", [self.model.homeostatic_objective[mol] for mol in fba.getHomeostaticTargetMolecules()])
		self.writeToListener("FBAResults", "homeostaticObjectiveValues", fba.getHomeostaticObjectiveValues())
		self.writeToListener("FBAResults", "kineticObjectiveValues", fba.getKineticObjectiveValues())
		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit", metabolite_counts_init)
		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal", metabolite_counts_final)
		self.writeToListener("EnzymeKinetics", "enzymeCountsInit", kinetic_enzyme_counts)
		self.writeToListener("EnzymeKinetics", "countsToMolar", counts_to_molar.asNumber(CONC_UNITS))
		self.writeToListener("EnzymeKinetics", "actualFluxes", fba.getReactionFluxes(self.model.kinetics_constrained_reactions) / time_step_unitless)
		self.writeToListener("EnzymeKinetics", "targetFluxes", targets / time_step_unitless)
		self.writeToListener("EnzymeKinetics", "targetFluxesUpper", upper_targets / time_step_unitless)
		self.writeToListener("EnzymeKinetics", "targetFluxesLower", lower_targets / time_step_unitless)
		self.writeToListener("EnzymeKinetics", "targetAAConc", [self.aa_targets.get(id_, 0) for id_ in self.aa_names])

	def update_amino_acid_targets(self, counts_to_molar):
		"""
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
		"""

		count_diff = self._sim.processes['PolypeptideElongation'].aa_count_diff

		if len(self.aa_targets):
			for aa, diff in count_diff.items():
				if aa in self.aa_targets_not_updated:
					continue
				self.aa_targets[aa] += diff
				# TODO (Santiago): Improve targets update
				if self.aa_targets[aa] < 0:
					print(f'Warning: updated amino acid target for {aa} was negative - adjusted to be positive.')
					self.aa_targets[aa] = 1

		# First time step of a simulation so set target to current counts to prevent
		# concentration jumps between generations
		else:
			for aa, counts in zip(self.aa_names, self.aas.total_counts()):
				if aa in self.aa_targets_not_updated:
					continue
				self.aa_targets[aa] = counts

		conc_updates = {aa: counts * counts_to_molar for aa, counts in self.aa_targets.items()}

		# Update linked metabolites that will follow an amino acid
		for met, link in self.linked_metabolites.items():
			conc_updates[met] = conc_updates.get(link['lead'], 0 * counts_to_molar) * link['ratio']

		return conc_updates

class FluxBalanceAnalysisModel(object):
	"""
	Metabolism model that solves an FBA problem with modular_fba.
	"""

	def __init__(self, sim_data, imports=None, timeline=None, include_ppgpp=True):
		# type: (SimulationDataEcoli, Optional[Set[str]], Optional[List[Tuple[float, str]]], bool) -> None
		"""
		Args:
			sim_data: simulation data
			imports: molecule IDs with compartment tag of molecules that can be
				imported into the cell
			timeline: timeline for nutrient changes during simulation
				(time of change, media ID), if None, nutrients for the saved
				condition are set at time 0 (eg. [(0.0, 'minimal')])
			include_ppgpp: if True, ppGpp is included as a concentration target
		"""

		if timeline is None:
			nutrients = sim_data.conditions[sim_data.condition]['nutrients']
			timeline = [(0., nutrients)]
		else:
			nutrients = timeline[0][1]

		if imports is None:
			imports = sim_data.external_state.exchange_data_from_media(nutrients)['importExchangeMolecules']

		# Local sim_data references
		metabolism = sim_data.process.metabolism
		constants = sim_data.constants
		mass = sim_data.mass

		# Load constants
		self.ngam = constants.non_growth_associated_maintenance
		gam = constants.darkATP * mass.cell_dry_mass_fraction

		self.exchange_constraints = metabolism.exchange_constraints

		self._biomass_concentrations = {}  # type: dict
		self._getBiomassAsConcentrations = mass.getBiomassAsConcentrations

		# Include ppGpp concentration target in objective if not handled kinetically in other processes
		self.ppgpp_id = sim_data.molecule_ids.ppGpp
		self.getppGppConc = sim_data.growth_rate_parameters.get_ppGpp_conc

		# go through all media in the timeline and add to metaboliteNames
		metaboliteNamesFromNutrients = set()
		if include_ppgpp:
			metaboliteNamesFromNutrients.add(self.ppgpp_id)
		for time, media_id in timeline:
			exchanges = sim_data.external_state.exchange_data_from_media(media_id)
			metaboliteNamesFromNutrients.update(
				metabolism.concentration_updates.concentrations_based_on_nutrients(
					imports=exchanges['importExchangeMolecules'])
				)
		self.metaboliteNamesFromNutrients = list(sorted(metaboliteNamesFromNutrients))
		exchange_molecules = sorted(sim_data.external_state.all_external_exchange_molecules)
		molecule_masses = dict(zip(exchange_molecules,
			sim_data.getter.get_masses(exchange_molecules).asNumber(MASS_UNITS / COUNTS_UNITS)))

		# Setup homeostatic objective concentration targets
		## Determine concentrations based on starting environment
		conc_dict = metabolism.concentration_updates.concentrations_based_on_nutrients(media_id=nutrients, imports=imports)
		doubling_time = sim_data.condition_to_doubling_time[sim_data.condition]
		conc_dict.update(self.getBiomassAsConcentrations(doubling_time))
		if include_ppgpp:
			conc_dict[self.ppgpp_id] = self.getppGppConc(doubling_time)
		self.homeostatic_objective = dict((key, conc_dict[key].asNumber(CONC_UNITS)) for key in conc_dict)

		## Include all concentrations that will be present in a sim for constant length listeners
		for met in self.metaboliteNamesFromNutrients:
			if met not in self.homeostatic_objective:
				self.homeostatic_objective[met] = 0.

		# Data structures to compute reaction bounds based on enzyme presence/absence
		self.catalyst_ids = metabolism.catalyst_ids
		self.reactions_with_catalyst = metabolism.reactions_with_catalyst

		i = metabolism.catalysis_matrix_I
		j = metabolism.catalysis_matrix_J
		v = metabolism.catalysis_matrix_V
		shape = (i.max() + 1, j.max() + 1)
		self.catalysis_matrix = csr_matrix((v, (i, j)), shape=shape)

		# Function to compute reaction targets based on kinetic parameters and molecule concentrations
		self.get_kinetic_constraints = metabolism.get_kinetic_constraints

		# Remove disabled reactions so they don't get included in the FBA problem setup
		kinetic_constraint_reactions = metabolism.kinetic_constraint_reactions
		constraintsToDisable = metabolism.constraints_to_disable
		self.active_constraints_mask = np.array([(rxn not in constraintsToDisable) for rxn in kinetic_constraint_reactions])
		self.kinetics_constrained_reactions = list(np.array(kinetic_constraint_reactions)[self.active_constraints_mask])

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
		fba_options = {
			"reactionStoich": metabolism.reaction_stoich,
			"externalExchangedMolecules": exchange_molecules,
			"objective": self.homeostatic_objective,
			"objectiveType": objective_type,
			"objectiveParameters": {
					"kineticObjectiveWeight": kinetic_objective_weight,
					'kinetic_objective_weight_in_range': kinetic_objective_weight_in_range,
					"reactionRateTargets": {reaction: 1 for reaction in self.kinetics_constrained_reactions},
					"oneSidedReactionTargets": [],
					},
			"moleculeMasses": molecule_masses,
			"secretionPenaltyCoeff": metabolism.secretion_penalty_coeff, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver": solver,
			"maintenanceCostGAM": gam.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction": metabolism.maintenance_reaction,
		}
		self.fba = FluxBalanceAnalysis(**fba_options)

		self.metabolite_names = {met: i for i, met in enumerate(self.fba.getOutputMoleculeIDs())}
		self.aa_names_no_location = [x[:-3] for x in sorted(sim_data.amino_acid_code_to_id_ordered.values())]

	def getBiomassAsConcentrations(self, doubling_time):
		"""
		Caches the result of the sim_data function to improve performance since
		function requires computation but won't change for a given doubling_time.

		Args:
			doubling_time (float with time units): doubling time of the cell to
				get the metabolite concentrations for

		Returns:
			dict {str : float with concentration units}: dictionary with metabolite
				IDs as keys and concentrations as values
		"""

		minutes = doubling_time.asNumber(units.min)  # hashable
		if minutes not in self._biomass_concentrations:
			self._biomass_concentrations[minutes] = self._getBiomassAsConcentrations(doubling_time)

		return self._biomass_concentrations[minutes]

	def update_external_molecule_levels(self, objective,
			metabolite_concentrations, external_molecule_levels):
		"""
		Limit amino acid uptake to what is needed to meet concentration objective
		to prevent use as carbon source, otherwise could be used as an infinite
		nutrient source.

		Args:
			objective (Dict[str, Unum]): homeostatic objective for internal
				molecules (molecule ID: concentration in counts/volume units)
			metabolite_concentrations (Unum[float]): concentration for each molecule
				in metabolite_names
			external_molecule_levels (np.ndarray[float]): current limits on external
				molecule availability

		Returns:
			external_molecule_levels (np.ndarray[float]): updated limits on external
				molcule availability

		TODO:
			- determine rate of uptake so that some amino acid uptake can
			be used as a carbon/nitrogen source
		"""

		external_exchange_molecule_ids = self.fba.getExternalMoleculeIDs()
		for aa in self.aa_names_no_location:
			if aa + "[p]" in external_exchange_molecule_ids:
				idx = external_exchange_molecule_ids.index(aa + "[p]")
			elif aa + "[c]" in external_exchange_molecule_ids:
				idx = external_exchange_molecule_ids.index(aa + "[c]")
			else:
				continue

			conc_diff = objective[aa + "[c]"] - metabolite_concentrations[self.metabolite_names[aa + "[c]"]].asNumber(CONC_UNITS)
			if conc_diff < 0:
				conc_diff = 0

			if external_molecule_levels[idx] > conc_diff:
				external_molecule_levels[idx] = conc_diff

		return external_molecule_levels

	def set_molecule_levels(self, metabolite_counts, counts_to_molar,
			coefficient, current_media_id, unconstrained, constrained, conc_updates, aa_uptake_package=None):
		"""
		Set internal and external molecule levels available to the FBA solver.

		Args:
			metabolite_counts (np.ndarray[int]): counts for each metabolite with a
				concentration target
			counts_to_molar (Unum): conversion from counts to molar (counts/volume units)
			coefficient (Unum): coefficient to convert from mmol/g DCW/hr to mM basis
				(mass.time/volume units)
			current_media_id (str): ID of current media
			unconstrained (Set[str]): molecules that have unconstrained import
			constrained (Dict[str, units.Unum]): molecules (keys) and their
				limited max uptake rates (values in mol / mass / time units)
			conc_updates (Dict[str, Unum]): updates to concentrations targets for
				molecules (molecule ID: concentration in mmol/L units)
			aa_uptake_package (Tuple[np.ndarray, np.ndarray, Boolean]): packed
				variables needed to set hard external amino acid updatkes
				(uptake rates, amino acid names, force levels)
		"""

		# Update objective from media exchanges
		external_molecule_levels, objective = self.exchange_constraints(
			self.fba.getExternalMoleculeIDs(), coefficient, CONC_UNITS,
			current_media_id, unconstrained, constrained, conc_updates,
			)
		self.fba.update_homeostatic_targets(objective)

		# Internal concentrations
		metabolite_conc = counts_to_molar * metabolite_counts
		self.fba.setInternalMoleculeLevels(metabolite_conc.asNumber(CONC_UNITS))

		# External concentrations
		external_molecule_levels = self.update_external_molecule_levels(
			objective, metabolite_conc, external_molecule_levels)
		self.fba.setExternalMoleculeLevels(external_molecule_levels)

		if aa_uptake_package:
			levels, molecules, force = aa_uptake_package
			self.fba.setExternalMoleculeLevels(levels, molecules=molecules, force=force)


	def set_reaction_bounds(self, catalyst_counts, counts_to_molar, coefficient,
			gtp_to_hydrolyze):
		"""
		Set reaction bounds for constrained reactions in the FBA object.

		Args:
			catalyst_counts (np.ndarray[int]): counts of enzyme catalysts
			counts_to_molar (Unum): conversion from counts to molar (counts/volume units)
			coefficient (Unum): coefficient to convert from mmol/g DCW/hr to mM basis
				(mass.time/volume units)
			gtp_to_hydrolyze (float): number of GTP molecules to hydrolyze to
				account for consumption in translation
		"""

		# Maintenance reactions
		## Calculate new NGAM
		flux = (self.ngam * coefficient).asNumber(CONC_UNITS)
		self.fba.setReactionFluxBounds(
			self.fba._reactionID_NGAM,
			lowerBounds=flux, upperBounds=flux,
			)

		## Calculate GTP usage based on how much was needed in polypeptide
		## elongation in previous step.
		flux = (counts_to_molar * gtp_to_hydrolyze).asNumber(CONC_UNITS)
		self.fba.setReactionFluxBounds(
			self.fba._reactionID_polypeptideElongationEnergy,
			lowerBounds=flux, upperBounds=flux,
			)

		# Set hard upper bounds constraints based on enzyme presence
		# (infinite upper bound) or absence (upper bound of zero)
		reaction_bounds = np.inf * np.ones(len(self.reactions_with_catalyst))
		no_rxn_mask = self.catalysis_matrix.dot(catalyst_counts) == 0
		reaction_bounds[no_rxn_mask] = 0
		self.fba.setReactionFluxBounds(self.reactions_with_catalyst,
			upperBounds=reaction_bounds, raiseForReversible=False)

	def set_reaction_targets(self, kinetic_enzyme_counts,
			kinetic_substrate_counts, counts_to_molar, time_step):
		# type: (np.ndarray, np.ndarray, units.Unum, units.Unum) -> Tuple[np.ndarray, np.ndarray, np.ndarray]
		"""
		Set reaction targets for constrained reactions in the FBA object.

		Args:
			kinetic_enzyme_counts (np.ndarray[int]): counts of enzymes used in
				kinetic constraints
			kinetic_substrate_counts (np.ndarray[int]): counts of substrates used
				in kinetic constraints
			counts_to_molar: conversion from counts to molar (float with counts/volume units)
			time_step: current time step (float with time units)

		Returns:
			mean_targets (np.ndarray[float]): mean target for each constrained reaction
			upper_targets (np.ndarray[float]): upper target limit for each constrained reaction
			lower_targets (np.ndarray[float]): lower target limit for each constrained reaction
		"""

		if self.use_kinetics:
			enzyme_conc = counts_to_molar * kinetic_enzyme_counts
			substrate_conc = counts_to_molar * kinetic_substrate_counts

			## Set target fluxes for reactions based on their most relaxed constraint
			reaction_targets = self.get_kinetic_constraints(enzyme_conc, substrate_conc)

			## Calculate reaction flux target for current time step
			targets = (time_step * reaction_targets).asNumber(CONC_UNITS)[
				self.active_constraints_mask, :]
			lower_targets = targets[:, 0]
			mean_targets = targets[:, 1]
			upper_targets = targets[:, 2]

			## Set kinetic targets only if kinetics is enabled
			self.fba.set_scaled_kinetic_objective(time_step.asNumber(units.s))
			self.fba.setKineticTarget(
				self.kinetics_constrained_reactions, mean_targets,
				lower_targets=lower_targets, upper_targets=upper_targets)
		else:
			lower_targets = np.zeros(len(self.kinetics_constrained_reactions))
			mean_targets = np.zeros(len(self.kinetics_constrained_reactions))
			upper_targets = np.zeros(len(self.kinetics_constrained_reactions))

		return mean_targets, upper_targets, lower_targets
