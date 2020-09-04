"""
SimulationData for metabolism process

TODO:
- improved estimate of ILE/LEU abundance or some external data point
- implement L1-norm minimization for AA concentrations
- find concentration for PI[c]
- add (d)NTP byproduct concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import absolute_import, division, print_function

from copy import copy
import re
from typing import Any, cast, Dict, Iterable, List, Optional, Set, Tuple, Union

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
# NOTE: Importing SimulationDataEcoli would make a circular reference so use Any.
#from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils import units
import six
from six.moves import range, zip


PPI_CONCENTRATION = 0.5e-3  # M, multiple sources
ECOLI_PH = 7.2
KINETIC_CONSTRAINT_CONC_UNITS = units.umol / units.L
K_CAT_UNITS = 1 / units.s
METABOLITE_CONCENTRATION_UNITS = units.mol / units.L

USE_ALL_CONSTRAINTS = False  # False will remove defined constraints from objective

REVERSE_TAG = ' (reverse)'
REVERSE_REACTION_ID = '{{}}{}'.format(REVERSE_TAG)
ENZYME_REACTION_ID = '{}__{}'

# Manually added concentration fold changes expected in media conditions
# relative to basal. Added in addition to the loaded flat file.
RELATIVE_CHANGES = {
	'minimal_minus_oxygen': {
		'CAMP[c]': 3.7,  # Unden, Duchenne. DOI: 10.1007/BF00415284 (Table 1, 2)
		},
	}

VERBOSE = False


class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._set_solver_values(sim_data.constants)
		self._buildBiomass(raw_data, sim_data)
		self._buildMetabolism(raw_data, sim_data)
		self._build_ppgpp_reactions(raw_data, sim_data)
		self._build_transport_reactions(raw_data, sim_data)

	def _set_solver_values(self, constants):
		"""
		Sets values to be used in the FBA solver.

		Attributes set:
			solver (str): solver ID, should match a value in modular_fba.py
			kinetic_objective_weight (float): weighting for the kinetic objective,
				1-weighting for the homeostatic objective
			kinetic_objective_weight_in_range (float): weighting for deviations
				from the kinetic target within min and max ranges
			secretion_penalty_coeff (float): penalty on secretion fluxes
		"""

		self.solver = "glpk-linear"
		if "linear" in self.solver:
			self.kinetic_objective_weight = constants.metabolismKineticObjectiveWeightLinear
		else:
			self.kinetic_objective_weight = constants.metabolismKineticObjectiveWeightQuadratic
		self.kinetic_objective_weight_in_range = constants.metabolism_kinetic_objective_weight_in_range
		self.secretion_penalty_coeff = constants.secretion_penalty_coeff

	def _buildBiomass(self, raw_data, sim_data):
		wildtypeIDs = set(entry["molecule id"] for entry in raw_data.biomass)

		# Create vector of metabolite target concentrations

		# Since the data only covers certain metabolites, we need to rationally
		# expand the dataset to include the other molecules in the biomass
		# function.

		# First, load in metabolites that do have concentrations, then assign
		# compartments according to those given in the biomass objective.  Or,
		# if there is no compartment, assign it to the cytoplasm.

		concentration_sources = [
			'Park Concentration',
			'Lempp Concentration',
			'Kochanowski Concentration',
			]
		excluded = {
			'Park Concentration': {
				'GLT',  # Steady state concentration reached with tRNA charging is much lower than Park
				},
			'Lempp Concentration': {
				'ATP',  # TF binding does not solve with average concentration
				},
			'Kochanowski Concentration': {
				'ATP',  # TF binding does not solve with average concentration
				},
			}
		metaboliteIDs = []
		metaboliteConcentrations = []

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			} # this assumes biomass reaction components only exist in a single compartment

		for row in raw_data.metaboliteConcentrations:
			metabolite_id = row['Metabolite']
			if not sim_data.getter.check_valid_molecule(metabolite_id):
				if VERBOSE:
					print('Metabolite concentration for unknown molecule: {}'
						.format(metabolite_id))
				continue

			# Use average of both sources
			# TODO (Travis): geometric mean?
			conc = np.nanmean([
				row[source].asNumber(METABOLITE_CONCENTRATION_UNITS)
				for source in concentration_sources
				if metabolite_id not in excluded.get(source, set())
				])

			# Check that a value was in the datasets being used
			if not np.isfinite(conc):
				if VERBOSE:
					print('No concentration in active datasets for {}'.format(metabolite_id))
				continue

			if metabolite_id in wildtypeIDtoCompartment:
				metaboliteIDs.append(
					metabolite_id + wildtypeIDtoCompartment[metabolite_id]
					)
			else:
				metaboliteIDs.append(
					metabolite_id + "[c]"
					)

			metaboliteConcentrations.append(conc)

		# CYS/SEL: concentration based on other amino acids
		aaConcentrations = []
		for aaIndex, aaID in enumerate(sim_data.amino_acid_code_to_id_ordered.values()):
			if aaID in metaboliteIDs:
				metIndex = metaboliteIDs.index(aaID)
				aaConcentrations.append(metaboliteConcentrations[metIndex])
		aaSmallestConc = min(aaConcentrations)

		metaboliteIDs.append("CYS[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		metaboliteIDs.append("L-SELENOCYSTEINE[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		# DGTP: set to smallest of all other DNTP concentrations
		dntpConcentrations = []
		for dntpIndex, dntpID in enumerate(sim_data.moleculeGroups.dntps):
			if dntpID in metaboliteIDs:
				metIndex = metaboliteIDs.index(dntpID)
				dntpConcentrations.append(metaboliteConcentrations[metIndex])
		dntpSmallestConc = min(dntpConcentrations)

		metaboliteIDs.append("DGTP[c]")
		metaboliteConcentrations.append(dntpSmallestConc)

		# H: from reported pH
		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append(sim_data.moleculeIds.proton)
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI: multiple sources report 0.5 mM
		metaboliteIDs.append(sim_data.moleculeIds.ppi)
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		metaboliteIDs.append("PI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		# include metabolites that are part of biomass
		for key, value in six.viewitems(sim_data.mass.getBiomassAsConcentrations(sim_data.doubling_time)):
			metaboliteIDs.append(key)
			metaboliteConcentrations.append(value.asNumber(METABOLITE_CONCENTRATION_UNITS))

		# Load relative metabolite changes
		relative_changes = {}
		for row in raw_data.relative_metabolite_concentrations:
			met = row['Metabolite']
			met_id = met + wildtypeIDtoCompartment.get(met, '[c]')

			# AA concentrations are determined through charging
			if met_id in sim_data.moleculeGroups.amino_acids:
				continue

			# Get relative metabolite change in each media condition
			for col, value in row.items():
				# Skip the ID column and minimal column (only has values of 1)
				# or skip invalid values
				if col == 'Metabolite' or col == 'minimal' or not np.isfinite(value):
					continue

				if col not in relative_changes:
					relative_changes[col] = {}
				relative_changes[col][met_id] = value

		## Add manually curated values for other media
		for media, data in RELATIVE_CHANGES.items():
			if media not in relative_changes:
				relative_changes[media] = {}
			for met, change in data.items():
				if met not in relative_changes[media]:
					relative_changes[media][met] = change

		# save concentrations as class variables
		unique_ids, counts = np.unique(metaboliteIDs, return_counts=True)
		if np.any(counts > 1):
			raise ValueError('Multiple concentrations for metabolite(s): {}'.format(', '.join(unique_ids[counts > 1])))

		# TODO (Travis): only pass raw_data and sim_data and create functions to load absolute and relative concentrations
		self.concentrationUpdates = ConcentrationUpdates(dict(zip(
			metaboliteIDs,
			METABOLITE_CONCENTRATION_UNITS * np.array(metaboliteConcentrations)
			)),
			relative_changes,
			raw_data.equilibriumReactions,
			sim_data.external_state.exchange_dict,
		)
		self.concDict = self.concentrationUpdates.concentrationsBasedOnNutrients("minimal")
		self.nutrientsToInternalConc = {}
		self.nutrientsToInternalConc["minimal"] = self.concDict.copy()

	def _buildMetabolism(self, raw_data, sim_data):
		"""
		Build the matrices/vectors for metabolism (FBA)
		Reads in and stores reaction and kinetic constraint information
		"""

		(reactionStoich, reversibleReactions, catalysts
			) = self.extract_reactions(raw_data, sim_data)

		# Load kinetic reaction constraints from raw_data
		known_metabolites = set(self.concDict)
		raw_constraints = self.extract_kinetic_constraints(raw_data, sim_data,
			stoich=reactionStoich, catalysts=catalysts,
			known_metabolites=known_metabolites)

		# Make modifications from kinetics data
		(constraints, reactionStoich, catalysts, reversibleReactions
			) = self._replace_enzyme_reactions(
			raw_constraints, reactionStoich, catalysts, reversibleReactions)

		# Create symbolic kinetic equations
		(self.kinetic_constraint_reactions, self.kinetic_constraint_enzymes,
			self.kinetic_constraint_substrates, self._kcats, self._saturations,
			self._enzymes, self.constraint_is_kcat_only
			) = self._lambdify_constraints(constraints)
		self._compiled_enzymes = None
		self._compiled_saturation = None

		# TODO: move this to a sim_data analysis script
		if VERBOSE:
			print('\nSummary of included metabolism kinetics:')
			print('Reactions with kinetics: {}'.format(len(self.kinetic_constraint_reactions)))
			print('Enzymes with kinetics: {}'.format(len(self.kinetic_constraint_enzymes)))
			print('Metabolites in kinetics: {}'.format(len(self.kinetic_constraint_substrates)))
			print('Number of kcat values: {}'.format(len([k for c in constraints.values() for k in c['kcat']])))
			print('Number of saturation terms: {}'.format(len([s for c in constraints.values() for s in c['saturation']])))

		# Verify no substrates with unknown concentrations have been added
		unknown = {m for m in self.kinetic_constraint_substrates
			if m not in known_metabolites}
		if unknown:
			raise ValueError('Unknown concentration for {}. Need to remove'
				' kinetics saturation term.'.format(', '.join(unknown)))

		# Extract data
		reactions_with_catalyst = sorted(catalysts)
		catalyst_ids = sorted({c for all_cat in catalysts.values()
			for c in all_cat})

		# Create catalysis matrix (to be used in the simulation)
		catalysisMatrixI = []
		catalysisMatrixJ = []
		catalysisMatrixV = []

		for row, reaction in enumerate(reactions_with_catalyst):
			for catalyst in catalysts[reaction]:
				col = catalyst_ids.index(catalyst)
				catalysisMatrixI.append(row)
				catalysisMatrixJ.append(col)
				catalysisMatrixV.append(1)

		catalysisMatrixI = np.array(catalysisMatrixI)
		catalysisMatrixJ = np.array(catalysisMatrixJ)
		catalysisMatrixV = np.array(catalysisMatrixV)

		# Properties for FBA reconstruction
		self.reactionStoich = reactionStoich
		self.maintenanceReaction = {"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "PI[c]": +1, "PROTON[c]": +1,}

		# Properties for catalysis matrix (to set hard bounds)
		self.reactionCatalysts = catalysts
		self.catalyst_ids = catalyst_ids
		self.reactions_with_catalyst = reactions_with_catalyst
		self.catalysisMatrixI = catalysisMatrixI
		self.catalysisMatrixJ = catalysisMatrixJ
		self.catalysisMatrixV = catalysisMatrixV

		# Properties for setting flux targets
		self.useAllConstraints = USE_ALL_CONSTRAINTS
		self.constraintsToDisable = [rxn["disabled reaction"]
			for rxn in raw_data.disabledKineticReactions]

	def _build_ppgpp_reactions(self, raw_data, sim_data):
		'''
		Creates structures for ppGpp reactions for use in polypeptide_elongation.

		Attributes set:
			ppgpp_synthesis_reaction (str): reaction ID for ppGpp synthesis
				(catalyzed by RelA and SpoT)
			ppgpp_degradation_reaction (str): reaction ID for ppGpp degradation
				(catalyzed by SpoT)
			ppgpp_reaction_names (list[str]): names of reaction involved in ppGpp
			ppgpp_reaction_metabolites (list[str]): names of metabolites in
				ppGpp reactions
			ppgpp_reaction_stoich (array[int]): 2D array with metabolites on rows
				and reactions on columns containing the stoichiometric coefficient
		'''

		self.ppgpp_synthesis_reaction = 'GDPPYPHOSKIN-RXN'
		self.ppgpp_degradation_reaction = 'PPGPPSYN-RXN'

		self.ppgpp_reaction_names = [
			self.ppgpp_synthesis_reaction,
			self.ppgpp_degradation_reaction,
			]

		self.ppgpp_reaction_metabolites = []

		# Indices (i: metabolite, j: reaction) and values (v: stoichiometry)
		# for sparse reaction matrix
		metabolite_indices = {}
		new_index = 0
		rxn_i = []
		rxn_j = []
		rxn_v = []

		# Record sparse indices in the matrix
		for j, rxn in enumerate(self.ppgpp_reaction_names):
			for met, stoich in self.reactionStoich[rxn].items():
				idx = metabolite_indices.get(met, new_index)

				if idx == new_index:
					metabolite_indices[met] = new_index
					self.ppgpp_reaction_metabolites.append(met)
					new_index += 1

				rxn_i.append(idx)
				rxn_j.append(j)
				rxn_v.append(stoich)

		# Assemble matrix based on indices
		# new_index is number of metabolites, j+1 is number of reactions
		self.ppgpp_reaction_stoich = np.zeros((new_index, j+1), dtype=np.int32)
		self.ppgpp_reaction_stoich[rxn_i, rxn_j] = rxn_v

	def _build_transport_reactions(self, raw_data, sim_data):
		"""
		Creates list of transport reactions that are included in the
		reaction network.

		Attributes set:
			transport_reactions (List[str]): transport reaction IDs in the
				metabolic network (includes reverse reactions and reactions
				with kinetic constraints)
		"""

		rxn_mapping = {}
		for rxn in sorted(self.reactionStoich):
			basename = rxn.split('__')[0].split(' (reverse)')[0]
			rxn_mapping[basename] = rxn_mapping.get(basename, []) + [rxn]

		self.transport_reactions = [
			rxn
			for row in raw_data.transport_reactions
			for rxn in rxn_mapping.get(row['reaction id'], [])
			]

	def getKineticConstraints(self, enzymes, substrates):
		# type: (units.Unum, units.Unum) -> units.Unum
		'''
		Allows for dynamic code generation for kinetic constraint calculation
		for use in Metabolism process. Inputs should be unitless but the order
		of magnitude should match the kinetics parameters (umol/L/s).

		If trying to pickle sim_data object after function has been called,
		_compiled_enzymes and _compiled_saturation might not be able to be pickled.
		See __getstate__(), __setstate__() comments on PR 111 to address.

		Returns np.array of floats of the kinetic constraint target for each
		reaction with kinetic parameters

		Args:
			enzymes: concentrations of enzymes associated with kinetic
				constraints (mol / volume units)
			substrates: concentrations of substrates associated with kinetic
				constraints (mol / volume units)

		Returns:
			(n reactions, 3): min, mean and max kinetic constraints for each
				reaction with kinetic constraints (mol / volume / time units)
		'''

		if self._compiled_enzymes is None:
			self._compiled_enzymes = eval('lambda e: {}'.format(self._enzymes))
		if self._compiled_saturation is None:
			self._compiled_saturation = eval('lambda s: {}'.format(self._saturations))

		# Strip units from args
		enzs = enzymes.asNumber(KINETIC_CONSTRAINT_CONC_UNITS)
		subs = substrates.asNumber(KINETIC_CONSTRAINT_CONC_UNITS)

		capacity = np.array(self._compiled_enzymes(enzs))[:, None] * self._kcats
		saturation = np.array([
			[min(v), sum(v) / len(v), max(v)]
			for v in self._compiled_saturation(subs)
			])

		return KINETIC_CONSTRAINT_CONC_UNITS * K_CAT_UNITS * capacity * saturation

	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits, currentNutrients, unconstrained, constrained, concModificationsBasedOnCondition = None):
		"""
		Called during Metabolism process
		Returns the homeostatic objective concentrations based on the current nutrients
		Returns levels for external molecules available to exchange based on the current nutrients
		"""

		newObjective = self.concentrationUpdates.concentrationsBasedOnNutrients(
			media_id=currentNutrients, conversion_units=targetUnits)
		if concModificationsBasedOnCondition is not None:
			newObjective.update(concModificationsBasedOnCondition)

		externalMoleculeLevels = np.zeros(len(exchangeIDs), np.float64)

		for index, moleculeID in enumerate(exchangeIDs):
			if moleculeID in unconstrained:
				externalMoleculeLevels[index] = np.inf
			elif moleculeID in constrained:
				externalMoleculeLevels[index] = (
					constrained[moleculeID] * coefficient
					).asNumber(targetUnits)
			else:
				externalMoleculeLevels[index] = 0.

		return externalMoleculeLevels, newObjective

	def set_supply_constants(self, sim_data):
		"""
		Sets constants to determine amino acid supply during translation.

		Args:
			sim_data (SimulationData object)

		Sets class attributes:
			KI_aa_synthesis (ndarray[float]): KI for each AA for synthesis
				portion of supply (in units of METABOLITE_CONCENTRATION_UNITS)
			KM_aa_export (ndarray[float]): KM for each AA for export portion
				of supply (in units of METABOLITE_CONCENTRATION_UNITS)
			fraction_supply_rate (float): fraction of AA supply that comes from
				a base synthesis rate
			fraction_import_rate (ndarray[float]): fraction of AA supply that
				comes from AA import if nutrients are present
			base_aa_conc (ndarray[float]): expected AA conc in basal condition
				(in units of METABOLITE_CONCENTRATION_UNITS)

		Assumptions:
			- Each internal amino acid concentration in 'minimal_plus_amino_acids'
			media is not lower than in 'minimal' media

		TODO (Travis):
			Base on measured KI and KM values.
			Add impact from synthesis enzymes and transporters.
			Better handling of concentration assumption
		"""

		aa_ids = sim_data.moleculeGroups.amino_acids
		conc = self.concentrationUpdates.concentrationsBasedOnNutrients

		aa_conc_basal = np.array([
			conc('minimal')[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
			for aa in aa_ids])
		aa_conc_aa_media = np.array([
			conc('minimal_plus_amino_acids')[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
			for aa in aa_ids])

		# Lower concentrations might produce strange rates (excess supply or
		# negative import when present externally) and constants so raise
		# to double check the implementation
		if not np.all(aa_conc_basal <= aa_conc_aa_media):
			aas = np.array(aa_ids)[np.where(aa_conc_basal > aa_conc_aa_media)]
			raise ValueError('Check that amino acid concentrations should be lower in amino acid media for {}'.format(aas))

		f_inhibited = sim_data.constants.fraction_supply_inhibited
		f_exported = sim_data.constants.fraction_supply_exported

		# Assumed units of METABOLITE_CONCENTRATION_UNITS for KI and KM
		self.KI_aa_synthesis = f_inhibited * aa_conc_basal / (1 - f_inhibited)
		self.KM_aa_export = (1 / f_exported - 1) * aa_conc_aa_media
		self.fraction_supply_rate = 1 - f_inhibited + aa_conc_basal / (self.KM_aa_export + aa_conc_basal)
		self.fraction_import_rate = 1 - (self.fraction_supply_rate + 1 / (1 + aa_conc_aa_media / self.KI_aa_synthesis) - f_exported)

	def aa_supply_scaling(self, aa_conc, aa_present):
		"""
		Called during polypeptide_elongation process
		Determine amino acid supply rate scaling based on current amino acid
		concentrations.

		Args:
			aa_conc (ndarray[float] with mol / volume units): internal
				concentration for each amino acid
			aa_present (ndarray[bool]): whether each amino acid is in the
				external environment or not

		Returns:
			ndarray[float]: scaling for the supply of each amino acid with
				higher supply rate if >1, lower supply rate if <1
		"""

		aa_conc = aa_conc.asNumber(METABOLITE_CONCENTRATION_UNITS)

		aa_supply = self.fraction_supply_rate
		aa_import = aa_present * self.fraction_import_rate
		aa_synthesis = 1 / (1 + aa_conc / self.KI_aa_synthesis)
		aa_export = aa_conc / (self.KM_aa_export + aa_conc)
		supply_scaling = aa_supply + aa_import + aa_synthesis - aa_export

		return supply_scaling

	@staticmethod
	def extract_reactions(raw_data, sim_data):
		# type: (KnowledgeBaseEcoli, Any) -> Tuple[Dict[str, Dict[str, int]], List[str], Dict[str, List[str]]]
		"""
		Extracts reaction data from raw_data to build metabolism reaction
		network with stoichiometry, reversibility and enzyme catalysts.

		Args:
			raw_data: knowledge base data
			sim_data (SimulationDataEcoli): simulation data

		Returns:
			reaction_stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction
			reversible_reactions: reaction IDs for reactions that have a reverse
				complement, does not have reverse tag
			reaction_catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich
		"""

		# Initialize variables to store reaction information
		reaction_stoich = {}
		reversible_reactions = []
		reaction_catalysts = {}

		# Load and parse reaction information from raw_data
		for reaction in cast(Any, raw_data).reactions:
			reaction_id = reaction["reaction id"]
			stoich = reaction["stoichiometry"]
			reversible = reaction["is reversible"]

			if len(stoich) <= 1:
				raise Exception("Invalid biochemical reaction: {}, {}".format(reaction_id, stoich))

			reaction_stoich[reaction_id] = stoich

			catalysts_for_this_rxn = []
			for catalyst in reaction["catalyzed by"]:
				try:
					catalysts_with_loc = catalyst + sim_data.getter.get_location_tag(catalyst)
					catalysts_for_this_rxn.append(catalysts_with_loc)
				# If we don't have the catalyst in our reconstruction, drop it
				except KeyError:
					if VERBOSE:
						print('Skipping catalyst {} for {} since it is not in the model'
							.format(catalyst, reaction_id))

			if len(catalysts_for_this_rxn) > 0:
				reaction_catalysts[reaction_id] = catalysts_for_this_rxn

			# Add the reverse reaction
			if reversible:
				reverse_reaction_id = REVERSE_REACTION_ID.format(reaction_id)
				reaction_stoich[reverse_reaction_id] = {
					moleculeID:-stoichCoeff
					for moleculeID, stoichCoeff in six.viewitems(reaction_stoich[reaction_id])
					}

				reversible_reactions.append(reaction_id)
				if len(catalysts_for_this_rxn) > 0:
					reaction_catalysts[reverse_reaction_id] = list(reaction_catalysts[reaction_id])

		return reaction_stoich, reversible_reactions, reaction_catalysts

	@staticmethod
	def match_reaction(stoich, catalysts, rxn_to_match, enz, mets, direction=None):
		# type: (Dict[str, Dict[str, int]], Dict[str, List[str]], str, str, List[str], Optional[str]) -> List[str]
		"""
		Matches a given reaction (rxn_to_match) to reactions that exist in
		stoich given that enz is known to catalyze the reaction and mets are
		reactants in the reaction. Can perform a fuzzy reaction match since
		rxn_to_match just needs to be part of the actual reaction name to match
		specific instances of a reaction.
		(eg. rxn_to_match="ALCOHOL-DEHYDROG-GENERIC-RXN" can match
		"ALCOHOL-DEHYDROG-GENERIC-RXN-ETOH/NAD//ACETALD/NADH/PROTON.30.").

		Args:
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction
			catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts,
				likely a subset of reactions in stoich
			rxn_to_match: reaction ID from kinetics to match to existing reactions
			enz: enzyme ID with location tag
			mets: metabolite IDs with no location tag from kinetics
			direction: reaction directionality, 'forward' or 'reverse' or None

		Returns:
			rxn_matches: matched reaction IDs in stoich
		"""

		# Mapping to handle instances of metabolite classes in kinetics
		# Keys: specific molecules in kinetics file
		# Values: class of molecules in reactions file that contain the key
		class_mets = {
			'RED-THIOREDOXIN-MONOMER': 'Red-Thioredoxin',
			'RED-THIOREDOXIN2-MONOMER': 'Red-Thioredoxin',
			'RED-GLUTAREDOXIN': 'Red-Glutaredoxins',
			'GRXB-MONOMER': 'Red-Glutaredoxins',
			'GRXC-MONOMER': 'Red-Glutaredoxins',
			'OX-FLAVODOXIN1': 'Oxidized-flavodoxins',
			'OX-FLAVODOXIN2': 'Oxidized-flavodoxins',
		}

		# Match full reaction name from partial reaction in kinetics. Must
		# also match metabolites since there can be multiple reaction instances.
		match = False
		match_candidates = []
		if rxn_to_match in stoich:
			match_candidates.append(rxn_to_match)
		else:
			for long_rxn, long_mets in stoich.items():
				if rxn_to_match in long_rxn and not long_rxn.endswith(REVERSE_TAG):
					match = True
					stripped_enzs = {e[:-3] for e in catalysts.get(long_rxn, [])}
					stripped_mets = {m[:-3] for m in long_mets}
					if (np.all([class_mets.get(m, m) in stripped_mets for m in mets])
							and enz in stripped_enzs):
						match_candidates.append(long_rxn)

		if len(match_candidates) == 0:
			if VERBOSE:
				if match:
					print('Partial reaction match: {} {} {} {} {}'.format(
						rxn_to_match, enz, stripped_enzs, mets, stripped_mets))
				else:
					print('No reaction match: {}'.format(rxn_to_match))

		# Determine direction of kinetic reaction from annotation or
		# metabolite stoichiometry.
		rxn_matches = []
		for rxn in match_candidates:
			reverse_rxn = REVERSE_REACTION_ID.format(rxn)
			reverse_rxn_exists = reverse_rxn in stoich
			if direction:
				reverse = direction == 'reverse'
			else:
				s = {k[:-3]: v for k, v in stoich.get(rxn, {}).items()}
				direction_ = np.unique(np.sign([
					s.get(class_mets.get(m, m), 0) for m in mets]))
				if len(direction_) == 0 and not reverse_rxn_exists:
					reverse = False
				elif len(direction_) != 1 or direction_[0] == 0:
					if VERBOSE:
						print('Conflicting directionality: {} {} {}'.format(
							rxn, mets, direction_))
					continue
				else:
					reverse = direction_[0] > 0

			# Verify a reverse reaction exists in the model
			if reverse:
				if reverse_rxn_exists:
					rxn_matches.append(reverse_rxn)
					continue
				else:
					if VERBOSE:
						print('No reverse reaction: {} {}'.format(rxn, mets))
					continue

			rxn_matches.append(rxn)

		return sorted(rxn_matches)

	@staticmethod
	def temperature_adjusted_kcat(kcat, temp=''):
		# type: (units.Unum, Union[float, str]) -> np.ndarray
		"""
		Args:
			kcat: enzyme turnover number(s) (1 / time)
			temp: temperature of measurement, defaults to 25 if ''

		Returns:
			temperature adjusted kcat values, in units of 1/s
		"""

		if isinstance(temp, str):
			temp = 25
		return 2**((37. - temp) / 10.) * kcat.asNumber(K_CAT_UNITS)

	@staticmethod
	def _construct_default_saturation_equation(mets, kms, kis, known_mets):
		# type: (List[str], List[float], List[float], Iterable[str]) -> str
		"""
		Args:
			mets: metabolite IDs with location tag for KM and KI
				parameters ordered to match order of kms then kis
			kms: KM parameters associated with mets
			kis: KI parameters associated with mets
			known_mets: metabolite IDs with location tag with known
				concentrations

		Returns:
			saturation equation with metabolites to replace delimited
				by double quote (eg. "metabolite")
		"""

		# Check input dimensions
		n_params = len(kms) + len(kis)
		if n_params == 0:
			return '1'
		if n_params != len(mets):
			if VERBOSE:
				print('Saturation parameter mismatch: {} {} {}'
					.format(mets, kms, kis))
			return '1'

		terms = []
		# Add KM terms
		for m, k in zip(mets, kms):
			if m in known_mets:
				terms.append('1+{}/"{}"'.format(k, m))
			elif VERBOSE:
				print('Do not have concentration for {} with KM={}'.format(m, k))
		# Add KI terms
		for m, k in zip(mets[len(kms):], kis):
			if m in known_mets:
				terms.append('1+"{}"/{}'.format(m, k))
			elif VERBOSE:
				print('Do not have concentration for {} with KI={}'.format(m, k))

		# Enclose groupings if being multiplied together
		if len(terms) > 1:
			terms[0] = '(' + terms[0]
			terms[-1] += ')'
		elif len(terms) == 0:
			return '1'

		return '1/({})'.format(')*('.join(terms))

	@staticmethod
	def _extract_custom_constraint(constraint, reactant_tags, product_tags, known_mets):
		# type: (Dict[str, Any], Dict[str, str], Dict[str, str], Set[str]) -> Tuple[Optional[np.ndarray], List[str]]
		"""
		Args:
			constraint: values defining a kinetic constraint with key:
				'customRateEquation' (str): mathematical representation of
					rate, must contain 'kcat*E'
				'customParameterVariables' (Dict[str, str]): mapping of
					variable names in the rate equation to metabolite IDs
					without location tags, must contain key 'E' (enzyme)
				'customParameterConstants' (List[str]): constant strings
					in the rate equation that correspond to values, must
					contain 'kcat'
				'customParameterConstantValues' (List[float]): values for
					each of the constant strings
				'Temp' (float or ''): temperature of measurement
			reactant_tags: mapping of molecule IDs without a location tag to
				molecule IDs with a location tag for all reactants
			product_tags: mapping of molecule IDs without a location tag to
				molecule IDs with a location tag for all products
			known_mets: molecule IDs with a location tag for molecules with
				known concentrations

		Returns:
			kcats: temperature adjusted kcat value, in units of 1/s
			saturation: saturation equation with metabolites to replace
				delimited by double quote (eg. "metabolite")
		"""

		equation = constraint['customRateEquation']
		variables = constraint['customParameterVariables']
		constant_keys = constraint['customParameterConstants']
		constant_values = constraint['customParameterConstantValues']
		temp = constraint['Temp']

		# Need to have these in the constraint
		kcat_str = 'kcat'
		enzyme_str = 'E'
		capacity_str = '{}*{}'.format(kcat_str, enzyme_str)

		# Need to replace these symbols in equations
		symbol_sub = {
			'^': '**',
			}

		# Make sure kcat exists
		if kcat_str not in constant_keys:
			if VERBOSE:
				print('Missing {} in custom constants: {}'.format(
					kcat_str, constant_keys))
			return None, []

		custom_kcat = 1 / units.s * np.array([constant_values[constant_keys.index(kcat_str)]])
		kcats = Metabolism.temperature_adjusted_kcat(custom_kcat, temp)

		# Make sure equation can be parsed, otherwise just return kcat
		if enzyme_str not in variables:
			if VERBOSE:
				print('Missing enzyme key ({}) in custom variables: {}'.format(
					enzyme_str, variables))
			return kcats, []
		if not capacity_str in equation:
			if VERBOSE:
				print('Expected to find {} in custom equation: {}'.format(
					capacity_str, equation))
			return kcats, []
		if len(constant_keys) != len(constant_values):
			if VERBOSE:
				print('Mismatch between constants: {} {}'.format(
					constant_keys, constant_values))
			return kcats, []

		variables_with_tags = {
			k: reactant_tags.get(v, product_tags.get(v, None))
			for k, v in variables.items()
			if k != enzyme_str and (v in reactant_tags or v in product_tags)
		}

		# Substitute values into custom equations
		## Replace terms with known constant values or sim molecule IDs with concentrations
		custom_subs = {k: str(v) for k, v in zip(constant_keys, constant_values)}
		custom_subs.update({
			k: '"{}"'.format(v)
			for k, v in variables_with_tags.items()
			if v in known_mets
		})

		## Remove capacity to get only saturation term
		new_equation = equation.replace(capacity_str, '1')

		## Tokenize equation to terms and symbols
		parsed_variables = re.findall(r'\w*', new_equation)[:-1]  # Remove trailing empty match
		## Ensure valid input of known variables or a float term
		for v in parsed_variables:
			if not (v == '' or v in custom_subs):
				try:
					float(v)
				except ValueError:
					if VERBOSE:
						print('Unknown value encountered in custom equation {}: {}'
							.format(equation, v))
					return kcats, []
		parsed_symbols = re.findall(r'\W', new_equation)
		tokenized_equation = np.array(parsed_variables)
		symbol_idx_mask = tokenized_equation == ''

		## Verify tokenized equation matches original before replacements
		tokenized_equation[symbol_idx_mask] = parsed_symbols
		if ''.join(tokenized_equation) != new_equation:
			if VERBOSE:
				print('Error parsing custom equation: {}'.format(equation))
			return kcats, []

		## Perform replacement of symbols
		tokenized_equation[symbol_idx_mask] = [symbol_sub.get(s, s) for s in parsed_symbols]

		# Reconstruct saturation equation with replacements
		saturation = [''.join([custom_subs.get(token, token) for token in tokenized_equation])]

		return kcats, saturation

	@staticmethod
	def extract_kinetic_constraints(raw_data, sim_data, stoich=None,
			catalysts=None, known_metabolites=None):
		# type: (KnowledgeBaseEcoli, Any, Optional[Dict[str, Dict[str, int]]], Optional[Dict[str, List[str]]], Optional[Set[str]]) -> Dict[Tuple[str, str], Dict[str, List[Any]]]
		"""
		Load and parse kinetic constraint information from raw_data

		Args:
			raw_data: knowledge base data
			sim_data (SimulationDataEcoli): simulation data
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction, if None, data
				is loaded from raw_data and sim_data
			catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich, if None, data is loaded from
				raw_data and sim_data
			known_metabolites: metabolites with known concentrations

		Returns:
			constraints: valid kinetic constraints for each reaction/enzyme pair
				{(reaction ID, enzyme with location tag): {
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}
		"""

		# Load data for optional args if needed
		if stoich is None or catalysts is None:
			loaded_stoich, _, loaded_catalysts = Metabolism.extract_reactions(raw_data, sim_data)

			if stoich is None:
				stoich = loaded_stoich
			if catalysts is None:
				catalysts = loaded_catalysts

		known_metabolites_ = set() if known_metabolites is None else known_metabolites

		constraints = {}  # type: Dict[Tuple[str, str], Dict[str, list]]
		for constraint in cast(Any, raw_data).metabolism_kinetics:
			rxn = constraint['reactionID']
			enzyme = constraint['enzymeID']
			metabolites = constraint['substrateIDs']
			direction = constraint['direction']
			kms = list(constraint['kM'].asNumber(KINETIC_CONSTRAINT_CONC_UNITS))
			kis = list(constraint['kI'].asNumber(KINETIC_CONSTRAINT_CONC_UNITS))
			n_reactants = len(metabolites) - len(kis)
			matched_rxns = Metabolism.match_reaction(stoich, catalysts, rxn, enzyme,
				metabolites[:n_reactants], direction)

			for matched_rxn in matched_rxns:
				# Ensure enzyme catalyzes reaction in model
				enzymes_tag_conversion = {e[:-3]: e for e in catalysts.get(matched_rxn, [])}
				if enzyme not in enzymes_tag_conversion:
					if VERBOSE:
						print('{} does not catalyze {}'.format(enzyme, matched_rxn))
					continue

				# Update metabolites with a location tag from the reaction
				# First look in reactants but some products can inhibit
				reactant_tags = {k[:-3]: k for k, v in stoich[matched_rxn].items() if v < 0}
				product_tags = {k[:-3]: k for k, v in stoich[matched_rxn].items() if v > 0}
				mets_with_tag = [
					reactant_tags.get(met, product_tags.get(met, ''))
					for met in metabolites
					if met in reactant_tags or met in product_tags
				]
				if len(mets_with_tag) != len(metabolites):
					# Warn if verbose but no continue since we can still use kcat
					if VERBOSE:
						print('Could not match all metabolites: {} {}'.format(
							metabolites, mets_with_tag))

				# Extract kcat and saturation parameters
				if constraint['rateEquationType'] == 'custom':
					kcats, saturation = Metabolism._extract_custom_constraint(
						constraint, reactant_tags, product_tags, known_metabolites_)
					if kcats is None:
						continue
				else:
					kcats = Metabolism.temperature_adjusted_kcat(constraint['kcat'], constraint['Temp'])
					if len(kcats) > 1:
						if len(kcats) != len(kms) or len(kms) != len(mets_with_tag):
							if VERBOSE:
								print('Could not align kcats and kms: {} {} {} {}'.format(
									rxn, kcats, kms, mets_with_tag))
							continue

						saturation = [
							Metabolism._construct_default_saturation_equation(
								[m], [km], [], known_metabolites_)
							for m, km in zip(mets_with_tag, kms)
						]
					else:
						saturation = [
							Metabolism._construct_default_saturation_equation(
								mets_with_tag, kms, kis, known_metabolites_)
						]

					saturation = [s for s in saturation if s != '1']

				# Add new kcats and saturation terms for the enzymatic reaction
				key = (matched_rxn, enzymes_tag_conversion[enzyme])
				entries = constraints.get(key, {})
				entries['kcat'] = entries.get('kcat', []) + list(kcats)
				entries['saturation'] = entries.get('saturation', []) + saturation
				constraints[key] = entries

		return constraints

	@staticmethod
	def _replace_enzyme_reactions(constraints, stoich, rxn_catalysts, reversible_rxns):
		# type: (Dict[Tuple[str, str], Dict[str, List[Any]]], Dict[str, Dict[str, int]], Dict[str, List[str]], List[str]) -> Tuple[Dict[str, Any], Dict[str, Dict[str, int]], Dict[str, List[str]], List[str]]
		"""
		Modifies reaction IDs in data structures to duplicate reactions with
		kinetic constraints and multiple enzymes.

		Args:
			constraints: valid kinetic constraints for each reaction/enzyme pair
				{(reaction ID, enzyme with location tag): {
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction, if None, data
				is loaded from raw_data and sim_data
			rxn_catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich, if None, data is loaded from
				raw_data and sim_data
			reversible_rxns: reaction IDs for reactions that have a reverse
				complement, does not have reverse tag

		Returns:
			new_constraints: valid kinetic constraints for each reaction
				{reaction ID: {
					'enzyme': enzyme catalyst (str),
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction with updated
				reactions for enzyme catalyzed kinetic reactions
			rxn_catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich with updated
				reactions for enzyme catalyzed kinetic reactions
			reversible_rxns: reaction IDs for reactions that have a reverse
				complement with updated reactions for enzyme catalyzed kinetic
				reactions, does not have reverse tag
		"""

		new_constraints = {}

		n_catalysts = {rxn: len(catalysts) for rxn, catalysts in rxn_catalysts.items()}

		# Split out reactions that are kinetically constrained and that have
		# more than one enzyme that catalyzes the reaction
		for (rxn, enzyme), constraint in constraints.items():
			if n_catalysts[rxn] > 1:
				# Create new reaction name with enzyme appended to the end
				if rxn.endswith(REVERSE_TAG):
					new_rxn = REVERSE_REACTION_ID.format(ENZYME_REACTION_ID.format(rxn[:-len(REVERSE_TAG)], enzyme[:-3]))
				else:
					new_rxn = ENZYME_REACTION_ID.format(rxn, enzyme[:-3])

				# Add the new reaction to appropriate lists and dicts
				stoich[new_rxn] = copy(stoich[rxn])
				rxn_catalysts[new_rxn] = [enzyme]
				if rxn in reversible_rxns:
					reversible_rxns.append(new_rxn)

				# Remove enzyme from old reaction and remove old reaction if no
				# more enzyme catalysts
				rxn_catalysts[rxn].pop(rxn_catalysts[rxn].index(enzyme))
				if len(rxn_catalysts[rxn]) == 0:
					stoich.pop(rxn)
					rxn_catalysts.pop(rxn)
					if rxn in reversible_rxns:
						reversible_rxns.pop(reversible_rxns.index(rxn))
			else:
				new_rxn = rxn

			# noinspection PyTypeChecker
			new_constraints[new_rxn] = dict(constraints[(rxn, enzyme)], enzyme=enzyme)

		return new_constraints, stoich, rxn_catalysts, reversible_rxns

	@staticmethod
	def _lambdify_constraints(constraints):
		# type: (Dict[str, Any]) -> Tuple[List[str], List[str], List[str], np.ndarray, str, str, np.ndarray]
		"""
		Creates str representations of kinetic terms to be used to create
		kinetic constraints that are returned with getKineticConstraints().

		Args:
			constraints: valid kinetic constraints for each reaction
				{reaction ID: {
					'enzyme': enzyme catalyst (str),
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}

		Returns:
			rxns: sorted reaction IDs for reactions with a kinetic constraint
			enzymes: sorted enzyme IDs for enzymes that catalyze a kinetic reaction
			substrates: sorted substrate IDs for substrates that are needed
				for kinetic saturation terms
			all_kcats: (n rxns, 3) min, mean and max kcat value for each reaction
			all_saturations: sympy str representation of a list of saturation
				terms (eg. '[s[0] / (1 + s[0]), 2 / (2 + s[1])]')
			all_enzymes: sympy str representation of enzymes for each reaction
				(eg. '[e[0], e[2], e[1]]')
			constraint_is_kcat_only: True if reaction only has kcat values and
				no saturation terms
		"""

		# Ordered lists of constraint related IDs
		rxns = sorted(constraints)
		enzymes = sorted({c['enzyme'] for c in constraints.values()})
		substrates = sorted({
			match.strip('"')
			for c in constraints.values()
			for s in c['saturation']
			for match in re.findall('".+?"', s)
		})

		# Mapping to replace molecule IDs with generic list strings
		enzyme_sub = {e: 'e[{}]'.format(i) for i, e in enumerate(enzymes)}
		substrate_sub = {'"{}"'.format(s): 's[{}]'.format(i)
			for i, s in enumerate(substrates)}

		# Mapping to replace generic list strings with sympy variables.
		# Need separate mapping from above because sympy handles '[]' as indexing
		# so location tags are not parsed properly.
		enzyme_symbols = {'e': [sp.symbols('e[{}]'.format(i)) for i in range(len(enzymes))]}
		substrate_symbols = {'s': [sp.symbols('s[{}]'.format(i)) for i in range(len(substrates))]}

		# Values to return
		all_kcats = np.zeros((len(rxns), 3))
		all_saturations = []
		all_enzymes = []
		constraint_is_kcat_only = []

		# Extract out data from each constraint
		for i, rxn in enumerate(rxns):
			kcats = constraints[rxn]['kcat']
			saturation = constraints[rxn]['saturation']
			enzyme = constraints[rxn]['enzyme']

			# Parse saturation equations into sympy format
			# If no saturation data is known, assume open range from 0 to 1
			saturations = []
			for sat in saturation:
				if sat == '1':
					continue

				for token, replace in substrate_sub.items():
					sat = sat.replace(token, replace)
				saturations.append(parse_expr(sat, local_dict=substrate_symbols))
			if len(saturations) == 0:
				saturations = [0, 1]
				constraint_is_kcat_only.append(True)
			else:
				constraint_is_kcat_only.append(False)

			# Save values for this constraint
			all_kcats[i, :] = [np.min(kcats), np.mean(kcats), np.max(kcats)]
			all_saturations.append(saturations)
			all_enzymes.append(parse_expr(enzyme_sub[enzyme], local_dict=enzyme_symbols))

		# Convert to str to save as class attr to be executed
		all_saturations = str(all_saturations)
		all_enzymes = str(all_enzymes)
		constraint_is_kcat_only = np.array(constraint_is_kcat_only)

		return rxns, enzymes, substrates, all_kcats, all_saturations, all_enzymes, constraint_is_kcat_only


# Class used to update metabolite concentrations based on the current nutrient conditions
class ConcentrationUpdates(object):
	def __init__(self, concDict, relative_changes, equilibriumReactions, exchange_data_dict):
		self.units = units.getUnit(list(concDict.values())[0])
		self.defaultConcentrationsDict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)
		self.exchange_fluxes = self._exchange_flux_present(exchange_data_dict)
		self.relative_changes = relative_changes

		# factor of internal amino acid increase if amino acids present in nutrients
		self.moleculeScaleFactors = {
			"L-ALPHA-ALANINE[c]": 2.,
			"ARG[c]": 2.,
			"ASN[c]": 2.,
			"L-ASPARTATE[c]": 2.,
			"CYS[c]": 2.,
			"GLT[c]": 1.1,
			"GLN[c]": 2.,
			"GLY[c]": 2.,
			"HIS[c]": 2.,
			"ILE[c]": 2.,
			"LEU[c]": 2.,
			"LYS[c]": 2.,
			"MET[c]": 2.,
			"PHE[c]": 2.,
			"PRO[c]": 2.,
			"SER[c]": 2.,
			"THR[c]": 2.,
			"TRP[c]": 2.,
			"TYR[c]": 2.,
			"L-SELENOCYSTEINE[c]": 2.,
			"VAL[c]": 2.,
		}

		self.moleculeSetAmounts = self._addMoleculeAmounts(equilibriumReactions, self.defaultConcentrationsDict)

	# return adjustments to concDict based on nutrient conditions
	def concentrationsBasedOnNutrients(self, media_id=None, conversion_units=None):
		if conversion_units:
			conversion = self.units.asNumber(conversion_units)
		else:
			conversion = self.units

		concentrationsDict = self.defaultConcentrationsDict.copy()

		metaboliteTargetIds = sorted(concentrationsDict.keys())
		concentrations = conversion * np.array([concentrationsDict[k] for k in metaboliteTargetIds])
		concDict = dict(zip(metaboliteTargetIds, concentrations))

		if media_id is not None:
			# For faster conversions than .asNumber(conversion_units) for each setAmount
			if conversion_units:
				conversion_to_no_units = conversion_units.asUnit(self.units)

			# Adjust for measured concentration changes in different media
			if media_id in self.relative_changes:
				for mol_id, conc_change in self.relative_changes[media_id].items():
					if mol_id in concDict:
						concDict[mol_id]  *= conc_change

			# Adjust for concentration changes based on presence in media
			exchanges = self.exchange_fluxes[media_id]
			for moleculeName, setAmount in six.viewitems(self.moleculeSetAmounts):
				if ((moleculeName in exchanges and (moleculeName[:-3] + "[c]" not in self.moleculeScaleFactors or moleculeName == "L-SELENOCYSTEINE[c]"))
						or (moleculeName in self.moleculeScaleFactors and moleculeName[:-3] + "[p]" in exchanges)):
					if conversion_units:
						setAmount = (setAmount / conversion_to_no_units).asNumber()
					concDict[moleculeName] = setAmount

		return concDict

	def _exchange_flux_present(self, exchange_data):
		# type: (Dict[str, Any]) -> Dict[str, Set[str]]
		"""
		Caches the presence of exchanges in each media condition based on
		exchange_data to set concentrations in concentrationsBasedOnNutrients().

		Args:
			exchange_data: dictionary of exchange data for all media conditions with keys:
				importUnconstrainedExchangeMolecules (dict[str, set[str]]): for each media ID key,
					exchange molecules (with location tag) that do not have an upper bound on their flux
				importConstrainedExchangeMolecules (dict[str, dict[str, float with mol/mass/time units]]):
					for each media ID key, constrained molecules (with location tag)
					with upper bound flux constraints

		Returns:
			sets of molecules IDs (with location tags) that can be imported for each
			media ID
		"""

		exchange_fluxes = {}
		for media, env in exchange_data.items():
			exchange_fluxes[media] = {mol for mol, conc in env.items() if conc > 0}

		return exchange_fluxes

	def _addMoleculeAmounts(self, equilibriumReactions, concDict):
		moleculeSetAmounts = {}
		for reaction in equilibriumReactions:
			# We only want to do this for species with standard Michaelis-Menten kinetics initially
			if len(reaction["stoichiometry"]) != 3:
				continue

			moleculeName = [x["molecule"] for x in reaction["stoichiometry"] if x["type"] == "metabolite"][0]
			amountToSet = 1e-4
			moleculeSetAmounts[moleculeName + "[p]"] = amountToSet * self.units
			moleculeSetAmounts[moleculeName + "[c]"] = amountToSet * self.units

		for moleculeName, scaleFactor in six.viewitems(self.moleculeScaleFactors):
			moleculeSetAmounts[moleculeName] = scaleFactor * concDict[moleculeName] * self.units
		return moleculeSetAmounts
