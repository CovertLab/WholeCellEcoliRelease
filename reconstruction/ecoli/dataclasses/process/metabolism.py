"""
SimulationData for metabolism process

TODO:
- improved estimate of ILE/LEU abundance or some external data point
- implement L1-norm minimization for AA concentrations
- find concentration for PI[c]
- add (d)NTP byproduct concentrations
"""

from __future__ import absolute_import, division, print_function

from copy import copy
import itertools
import re
from typing import Any, cast, Dict, Iterable, List, Optional, Set, Tuple, Union

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from reconstruction.ecoli.dataclasses.getter_functions import UNDEFINED_COMPARTMENT_IDS_TO_ABBREVS
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
# NOTE: Importing SimulationDataEcoli would make a circular reference so use Any.
#from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils import units
import six
from six.moves import range, zip

KINETIC_CONSTRAINT_CONC_UNITS = units.umol / units.L
K_CAT_UNITS = 1 / units.s
METABOLITE_CONCENTRATION_UNITS = units.mol / units.L
DRY_MASS_UNITS = units.fg

USE_ALL_CONSTRAINTS = False  # False will remove defined constraints from objective

REVERSE_TAG = ' (reverse)'
REVERSE_REACTION_ID = '{{}}{}'.format(REVERSE_TAG)
ENZYME_REACTION_ID = '{}__{}'

VERBOSE = False


class InvalidReactionDirectionError(Exception):
	pass


class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._set_solver_values(sim_data.constants)
		self._build_biomass(raw_data, sim_data)
		self._build_metabolism(raw_data, sim_data)
		self._build_ppgpp_reactions(raw_data, sim_data)
		self._build_transport_reactions(raw_data, sim_data)
		self._build_amino_acid_pathways(raw_data, sim_data)

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
			self.kinetic_objective_weight = constants.metabolism_kinetic_objective_weight_linear
		else:
			self.kinetic_objective_weight = constants.metabolism_kinetic_objective_weight_quadratic
		self.kinetic_objective_weight_in_range = constants.metabolism_kinetic_objective_weight_in_range
		self.secretion_penalty_coeff = constants.secretion_penalty_coeff

	def _build_biomass(self, raw_data, sim_data):
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
			'Sander Concentration',
			]
		excluded = {
			'Park Concentration': {
				'GLT',  # Steady state concentration reached with tRNA charging is much lower than Park
				'THR',  # Attenuation needs concentration to be lower to match validation data
				'VAL',  # Synthesis pathway kcat needs concentration to be lower and closer to KI
				},
			'Lempp Concentration': {
				'ATP',  # TF binding does not solve with average concentration
				'VAL',  # Synthesis pathway kcat needs concentration to be lower and closer to KI
				},
			'Kochanowski Concentration': {
				'ATP',  # TF binding does not solve with average concentration
				},
			'Sander Concentration': {
				'GLT',  # Steady state concentration reached with tRNA charging is much lower than Sander
				},
			}
		metaboliteIDs = []
		metaboliteConcentrations = []

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			} # this assumes biomass reaction components only exist in a single compartment

		for row in raw_data.metabolite_concentrations:
			metabolite_id = row['Metabolite']
			if not sim_data.getter.is_valid_molecule(metabolite_id):
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
		for dntpIndex, dntpID in enumerate(sim_data.molecule_groups.dntps):
			if dntpID in metaboliteIDs:
				metIndex = metaboliteIDs.index(dntpID)
				dntpConcentrations.append(metaboliteConcentrations[metIndex])
		dntpSmallestConc = min(dntpConcentrations)

		metaboliteIDs.append("DGTP[c]")
		metaboliteConcentrations.append(dntpSmallestConc)

		# H: from reported pH
		hydrogenConcentration = 10**(-sim_data.constants.pH)

		metaboliteIDs.append(sim_data.molecule_ids.proton)
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI
		ppi_conc = sim_data.constants.ppi_concentration.asNumber(
			METABOLITE_CONCENTRATION_UNITS)
		metaboliteIDs.append(sim_data.molecule_ids.ppi)
		metaboliteConcentrations.append(ppi_conc)

		metaboliteIDs.append("Pi[c]")
		metaboliteConcentrations.append(ppi_conc)

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
			if met_id in sim_data.molecule_groups.amino_acids:
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
		for media, data in sim_data.adjustments.relative_metabolite_concentrations_changes.items():
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
		conc_dict = dict(zip(metaboliteIDs, METABOLITE_CONCENTRATION_UNITS * np.array(metaboliteConcentrations)))
		all_metabolite_ids = {met['id'] for met in raw_data.metabolites}
		linked_metabolites = self._build_linked_metabolites(raw_data, conc_dict)
		self.concentration_updates = ConcentrationUpdates(
			conc_dict,
			relative_changes,
			raw_data.equilibrium_reactions,
			sim_data.external_state.exchange_dict,
			all_metabolite_ids,
			linked_metabolites,
		)
		self.conc_dict = self.concentration_updates.concentrations_based_on_nutrients("minimal")
		self.nutrients_to_internal_conc = {}
		self.nutrients_to_internal_conc["minimal"] = self.conc_dict.copy()

	def _build_linked_metabolites(self, raw_data, conc_dict):
		"""
		Calculates ratio between linked metabolites to keep it constant
		throughout a simulation.

		Returns:
			linked_metabolites (Dict[str, Dict[str, Any]]): mapping from a
				linked metabolite to its lead metabolite and concentration
				ratio to be maintained with the following keys:
					'lead' (str): metabolite to link the concentration to
					'ratio' (float): ratio to multiply the lead concentration by
		"""

		linked_metabolites = {}
		for row in raw_data.linked_metabolites:
			lead = row['Lead metabolite']
			linked = row['Linked metabolite']
			ratio = units.strip_empty_units(conc_dict[linked] / conc_dict[lead])

			linked_metabolites[linked] = {'lead': lead, 'ratio': ratio}

		return linked_metabolites

	def _build_metabolism(self, raw_data, sim_data):
		"""
		Build the matrices/vectors for metabolism (FBA)
		Reads in and stores reaction and kinetic constraint information
		"""

		(reaction_stoich, reversible_reactions, catalysts
			) = self.extract_reactions(raw_data, sim_data)

		# Load kinetic reaction constraints from raw_data
		known_metabolites = set(self.conc_dict)
		raw_constraints = self.extract_kinetic_constraints(raw_data, sim_data,
			stoich=reaction_stoich, catalysts=catalysts,
			known_metabolites=known_metabolites)

		# Make modifications from kinetics data
		(constraints, reaction_stoich, catalysts, reversible_reactions
			) = self._replace_enzyme_reactions(
			raw_constraints, reaction_stoich, catalysts, reversible_reactions)

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
		self.reaction_stoich = reaction_stoich
		# TODO (ggsun): add this as a raw .tsv file
		self.maintenance_reaction = {"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "Pi[c]": +1, "PROTON[c]": +1, }

		# Properties for catalysis matrix (to set hard bounds)
		self.reaction_catalysts = catalysts
		self.catalyst_ids = catalyst_ids
		self.reactions_with_catalyst = reactions_with_catalyst
		self.catalysis_matrix_I = catalysisMatrixI
		self.catalysis_matrix_J = catalysisMatrixJ
		self.catalysis_matrix_V = catalysisMatrixV

		# Properties for setting flux targets
		self.use_all_constraints = USE_ALL_CONSTRAINTS
		self.constraints_to_disable = [rxn["disabled reaction"]
			for rxn in raw_data.disabled_kinetic_reactions]

		self.amino_acid_export_kms = raw_data.amino_acid_export_kms

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
			for met, stoich in self.reaction_stoich[rxn].items():
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
		transport_reactions = [
			rxn_id for rxn_id, stoich in self.reaction_stoich.items()
			if self._is_transport_rxn(stoich)]

		self.transport_reactions = transport_reactions

	def _build_amino_acid_pathways(self, raw_data, sim_data):
		"""
		Creates mapping between enzymes and amino acid pathways with
		allosteric inhibition feedback from the amino acid.

		Attributes set:
			aa_synthesis_pathways (Dict[str, Dict]): data for allosteric
				inhibition of amino acid pathways indexed by amino acid ID with
				location tag and nested dictionary with the following keys:
					'enzymes' (str): limiting/regulated enzyme ID in synthesis
						pathway with location tag
					'kcat_data' (units.Unum): kcat associated with enzyme
						reaction with units of 1/time
					'ki' (Tuple[units.Unum, units.Unum]]): lower and upper
						limits of KI associated with enzyme reaction with units
						of mol/volume
		"""

		self.aa_synthesis_pathways = {}
		cytoplasm_tag = '[c]'

		for row in raw_data.amino_acid_pathways:
			data = {}
			data['enzymes'] = [e + cytoplasm_tag for e in row['Enzymes']]
			data['kcat_data'] = 0 / units.s if units.isnan(row['kcat']) else row['kcat']
			if units.isnan(row['KI, lower bound']) or units.isnan(row['KI, lower bound']):
				data['ki'] = None
			else:
				data['ki'] = (row['KI, lower bound'], row['KI, upper bound'])
			data['upstream'] = {k + cytoplasm_tag: v for k, v in row['Upstream amino acids'].items()}
			data['reverse'] = {k + cytoplasm_tag: v for k, v in row['Reverse amino acids'].items()}
			data['km, upstream'] = {k + cytoplasm_tag: v for k, v in row['KM, upstream'].items()}
			data['km, reverse'] = row['KM, reverse']
			data['km, degradation'] = np.inf * units.mol/units.L if units.isnan(row['KM, degradation']) else row['KM, degradation']
			data['downstream'] = {k + cytoplasm_tag: v for k, v in row['Downstream amino acids'].items()}
			self.aa_synthesis_pathways[row['Amino acid'] + cytoplasm_tag] = data

		self.aa_synthesis_pathway_adjustments = {}
		for row in raw_data.adjustments.amino_acid_pathways:
			# Read data from row
			aa = row['Amino acid'] + cytoplasm_tag
			parameter = row['Parameter']
			factor = row['Factor']

			# Store adjustments to be used later
			adjustments = self.aa_synthesis_pathway_adjustments.get(aa, {})
			adjustments[parameter] = factor
			self.aa_synthesis_pathway_adjustments[aa] = adjustments


	def get_kinetic_constraints(self, enzymes, substrates):
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

	def exchange_constraints(self, exchangeIDs, coefficient, targetUnits, media_id,
			unconstrained, constrained, concModificationsBasedOnCondition = None):
		"""
		Called during Metabolism process
		Returns the homeostatic objective concentrations based on the current nutrients
		Returns levels for external molecules available to exchange based on the current nutrients
		"""

		newObjective = self.concentration_updates.concentrations_based_on_nutrients(
			imports=unconstrained.union(constrained), media_id=media_id,
			conversion_units=targetUnits)
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

	def set_phenomological_supply_constants(self, sim_data):
		"""
		Sets constants to determine amino acid supply during translation.  Used
		with aa_supply_scaling() during simulations but supply can
		alternatively be determined mechanistically.  This approach may require
		manually adjusting constants (fraction_supply_inhibited and
		fraction_supply_exported) but has less variability related to gene
		expression and regulation.

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

		Assumptions:
			- Each internal amino acid concentration in 'minimal_plus_amino_acids'
			media is not lower than in 'minimal' media

		TODO (Travis):
			Better handling of concentration assumption
		"""

		aa_ids = sim_data.molecule_groups.amino_acids
		conc = self.concentration_updates.concentrations_based_on_nutrients

		aa_conc_basal = np.array([
			conc(media_id='minimal')[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
			for aa in aa_ids])
		aa_conc_aa_media = np.array([
			conc(media_id='minimal_plus_amino_acids')[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
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

	def get_aa_to_transporters_mapping_data(self, sim_data, export=False):
		'''
		Creates a dictionary that maps amino acids with their transporters.
		Based on this dictionary, it creates a correlation matrix with rows
		as AA and columns as transporters.

		Args:
			sim_data (SimulationData object)
			export (Boolean): if True, the parameters calculated are for mechanistic
				export instead of uptake

		Returns:
			aa_to_transporters (Dict[str, list]): dictonary that maps aa to
				transporters involved in transport reactions
			aa_to_transporters_matrix (np.ndarray[int]): correlation matrix.
				Columns correspond to transporter enzymes and rows to amino acids
			aa_transporters_names (np.ndarray[str]): names of all transporters
		'''

		def matches_direction(direction):
			if export:
				return direction < 0
			else:
				return direction > 0


		# Mapping aminoacids to their transporters
		# CYS does not have any uptake reaction, so we initialize the dict with it to ensure
		# the presence of the 21 AAs
		# TODO (Santiago): Reversible reactions?
		aa_to_transporters = {"CYS[c]": []}
		for reaction in self.transport_reactions:
			for aa in sim_data.molecule_groups.amino_acids:
				if aa in self.reaction_stoich[reaction] and matches_direction(self.reaction_stoich[reaction][aa]):
					if aa not in aa_to_transporters:
						aa_to_transporters[aa] = []
					aa_to_transporters[aa] += self.reaction_catalysts[reaction]

		aa_to_transporters = {aa: aa_to_transporters[aa] for aa in sim_data.molecule_groups.amino_acids}

		c = 0
		transporters_to_idx = {}
		for aa, transporters in aa_to_transporters.items():
			for transporter in transporters:
				if transporter not in transporters_to_idx:
					transporters_to_idx[transporter]=c
					c += 1

		aa_to_transporters_matrix = [0]*len(aa_to_transporters)

		for i, trnspts in enumerate(aa_to_transporters.values()):
			temp = [0] * len(transporters_to_idx)
			for tr in trnspts:
				temp[transporters_to_idx[tr]] = 1
			aa_to_transporters_matrix[i] = temp

		aa_transporters_names = list(transporters_to_idx.keys())

		return aa_to_transporters, np.array(aa_to_transporters_matrix), np.array(aa_transporters_names)

	def set_mechanistic_export_constants(self, sim_data, cell_specs, basal_container):
		'''
		Calls get_aa_to_transporters_mapping_data() for AA export, which calculates
		the total amount of export transporter counts per AA. Kcats are calculated using
		the same exchange rates as for uptake and transporter counts. Missing KMs are calculated
		based on present KMs. This is done by calculating the average factor for
		KMs compared to estimated concentration (av_factor = sum(KM / concentration) / n_aa_with_kms).
		** KM = av_factor * concentration


		Args:
			sim_data (SimulationData object)
			cell_specs (Dict[str, Dict])
			basal_container (BulkObjectsContainer): average initialization
				container in the basal condition

		Sets class attribute:
			aa_to_export_transporters (Dict[str, list]): dictonary that maps aa to
				transporters involved in export reactions
			aa_to_export_transporters_matrix (np.ndarray[int]): correlation matrix.
				Columns correspond to exporting enzymes and rows to amino acids
			aa_export_transporters_names (np.ndarray[str]): names of all exporters
			aa_export_kms (np.ndarray[float]): kms corresponding to generic transport/enzyme
				reactions for each AA in concentration units (METABOLITE_CONCENTRATION_UNITS) [M]
			export_kcats_per_aa (np.ndarray[float]): kcats corresponding to generic export
				reactions for each AA. Units in counts/second [1/s]
		'''

		self.aa_to_export_transporters, self.aa_to_export_transporters_matrix, self.aa_export_transporters_names = self.get_aa_to_transporters_mapping_data(
			sim_data, export=True)

		aa_names = sim_data.molecule_groups.amino_acids
		counts_to_molar = (sim_data.constants.cell_density / cell_specs['basal']['avgCellDryMassInit']) / sim_data.constants.n_avogadro
		aa_conc = {aa: counts * counts_to_molar for aa, counts in zip(aa_names, basal_container.counts(aa_names))}
		aa_with_km = {}

		# Calculate average factor to estimate missing KMs
		coeff_estimate_kms = 0
		for export_kms in self.amino_acid_export_kms:
			aa_with_km[export_kms['Amino Acid']] = export_kms['KM']
			coef_per_aa = 0
			for km in export_kms['KM'].values():
				coef_per_aa += km.asUnit(METABOLITE_CONCENTRATION_UNITS) / aa_conc[export_kms['Amino Acid']]
			coeff_estimate_kms += coef_per_aa / len(export_kms['KM'])
		coeff_estimate_kms = coeff_estimate_kms / len(self.amino_acid_export_kms)

		# Calculate estimated KMs for each AA
		single_kms = {}
		for aa in aa_names:
			if aa in aa_with_km:
				single_kms[aa] = np.mean(list(aa_with_km[aa].values()))
			else:
				single_kms[aa] = coeff_estimate_kms * aa_conc[aa]

		self.aa_export_kms = np.array([single_kms[aa].asNumber(METABOLITE_CONCENTRATION_UNITS) for aa in aa_names])

	def set_mechanistic_uptake_constants(self, sim_data, cell_specs, with_aa_container):
		'''
		Based on the matrix calculated in get_aa_to_transporters_mapping_data(), we calculate
		the total amount of transporter counts per AA.

		Args:
			sim_data (SimulationData object)
			cell_specs (Dict[str, Dict])
			with_aa_container (BulkObjectsContainer): average initialization
				container in the with_aa condition

		Sets class attribute:
			aa_to_transporters (Dict[str, list]): dictonary that maps aa to
				transporters involved in import reactions
			aa_to_transporters_matrix (np.ndarray[int]): correlation matrix.
				Columns correspond to importing enzymes and rows to amino acids
			aa_transporters_names (np.ndarray[str]): names of all importers
			uptake_kcats_per_aa (np.ndarray[float]): kcats corresponding to generic uptake
				reactions for each AA. Units in counts/second [1/s]

		TODO:
			- Include external amino acid concentrations and KM values
		'''

		aa_names = sim_data.molecule_groups.amino_acids
		counts_to_molar = (sim_data.constants.cell_density / cell_specs['with_aa']['avgCellDryMassInit']) / sim_data.constants.n_avogadro
		aa_counts = with_aa_container.counts(aa_names)
		exchange_rates = sim_data.process.metabolism.specific_import_rates * cell_specs['with_aa']['avgCellDryMassInit'].asNumber(units.fg)

		self.aa_to_transporters, self.aa_to_transporters_matrix, self.aa_transporters_names = self.get_aa_to_transporters_mapping_data(sim_data)

		importer_counts = with_aa_container.counts(self.aa_transporters_names)
		exporter_counts = with_aa_container.counts(self.aa_export_transporters_names)
		counts_per_aa_import = self.aa_to_transporters_matrix.dot(importer_counts)
		counts_per_aa_export = self.aa_to_export_transporters_matrix.dot(exporter_counts)
		kms = self.aa_export_kms / counts_to_molar.asNumber(METABOLITE_CONCENTRATION_UNITS)

		# Calculate kcats based on specific_import_rates, dry mass, transporters counts, export kms and counts of aas
		with np.errstate(divide='ignore'):
			vmax = exchange_rates / (1 - (aa_counts/(kms + aa_counts)))
			self.uptake_kcats_per_aa = vmax / counts_per_aa_import
			self.export_kcats_per_aa = vmax / counts_per_aa_export
		self.export_kcats_per_aa[counts_per_aa_export == 0] = 0
		self.uptake_kcats_per_aa[counts_per_aa_import == 0] = 0

	def set_mechanistic_supply_constants(self, sim_data, cell_specs, basal_container, with_aa_container):
		"""
		Sets constants to determine amino acid supply during translation.  Used
		with amino_acid_synthesis() and amino_acid_import() during simulations
		but supply can alternatively be determined phenomologically.  This
		approach is more detailed and should better respond to environmental
		changes and perturbations but has more variability related to gene
		expression and regulation.

		Args:
			sim_data (SimulationData object)
			cell_specs (Dict[str, Dict])
			basal_container (BulkObjectsContainer): average initialization
				container in the basal condition
			with_aa_container (BulkObjectsContainer): average initialization
				container in the with_aa condition

		Sets class attributes:
			aa_enzymes (np.ndarray[str]): enzyme ID with location tag for each
				enzyme that can catalyze an amino acid pathway with
				self.enzyme_to_amino_acid mapping these to each amino acid
			aa_kcats (np.ndarray[float]): kcat value for each synthesis pathway
				in units of K_CAT_UNITS, ordered by amino acid molecule group
			aa_kis (np.ndarray[float]): KI value for each synthesis pathway
				in units of METABOLITE_CONCENTRATION_UNITS, ordered by amino
				acid molecule group. Will be inf if there is no inhibitory
				control.
			aa_upstream_kms (List[List[float]]): KM value associated with the
				amino acid that feeds into each synthesis pathway in units of
				METABOLITE_CONCENTRATION_UNITS, ordered by amino acid molecule
				group. Will be 0 if there is no upstream amino acid considered.
			aa_reverse_kms (np.ndarray[float]): KM value associated with the
				amino acid in each synthesis pathway in units of
				METABOLITE_CONCENTRATION_UNITS, ordered by amino acid molecule
				group. Will be inf if the synthesis pathway is not reversible.
			aa_upstream_mapping (np.ndarray[int]): index of the upstream amino
				acid that feeds into each synthesis pathway, ordered by amino
				acid molecule group
			enzyme_to_amino_acid (np.ndarray[float]): relationship mapping from
				aa_enzymes to amino acids (n enzymes, m amino acids).  Will
				contain a 1 if the enzyme associated with the row can catalyze
				the pathway for the amino acid associated with the column
			aa_forward_stoich (np.ndarray[float]): relationship mapping from
				upstream amino acids to downstream amino acids (n upstream,
				m downstream).  Will contain a -1 if the amino acid associated
				with the row is required for synthesis of the amino acid
				associated with the column
			aa_reverse_stoich (np.ndarray[float]): relationship mapping from
				upstream amino acids to downstream amino acids (n downstream,
				m upstream).  Will contain a -1 if the amino acid associated
				with the row is produced through a reverse reaction from
				the amino acid associated with the column
			specific_import_rates (np.ndarray[float]): import rates expected
				in rich media conditions for each amino acid normalized by dry
				cell mass in units of K_CAT_UNITS / DRY_MASS_UNITS,
				ordered by amino acid molecule group

		Assumptions:
			- Only one reaction is limiting in an amino acid pathway (typically
			the first and one with KI)
			- kcat applies to forward and reverse reaction and multiple enzymes
			catalyzing the same pathway will have the same kcat and saturation
			terms

		TODO:
			Handle different kcats (or enzymes) for reverse reactions
			Search for new kcat/KM values in literature or use metabolism_kinetics.tsv
			Consider multiple reaction steps
		"""

		aa_ids = sim_data.molecule_groups.amino_acids
		conc = self.concentration_updates.concentrations_based_on_nutrients

		# Allosteric inhibition constants to match required supply rate
		rates = (
			sim_data.translation_supply_rate['minimal']
			* sim_data.mass.avg_cell_dry_mass_init * sim_data.constants.n_avogadro
			)
		supply = {
			aa: rate
			for aa, rate in zip(sim_data.molecule_groups.amino_acids, rates)
			}
		aa_enzymes = []
		enzyme_to_aa = []
		aa_kcats = {}
		aa_kis = {}
		upstream_aas_for_km = {}
		upstream_aas = {}
		reverse_aas = {}
		aa_upstream_kms = {}
		aa_reverse_kms = {}
		aa_degradation_kms = {}
		degradation_rates = {}
		minimal_conc = conc('minimal')

		# Get order of amino acids to calculate parameters for to ensure that
		# parameters that are dependent on other amino acids are run after
		# those calculations have completed
		dependencies = {}
		for aa in aa_ids:
			for downstream_aa in self.aa_synthesis_pathways[aa]['downstream']:
				if np.isfinite(self.aa_synthesis_pathways[downstream_aa]['km, degradation'].asNumber()):
					dependencies[aa] = dependencies.get(aa, []) + [downstream_aa]
		ordered_aa_ids = []
		for _ in aa_ids:  # limit number of iterations number of amino acids in case there are cyclic links
			for aa in sorted(set(aa_ids) - set(ordered_aa_ids)):
				for downstream_aa in dependencies.get(aa, []):
					if downstream_aa not in ordered_aa_ids:
						break
				else:
					ordered_aa_ids.append(aa)
		if len(ordered_aa_ids) != len(aa_ids):
			raise RuntimeError('Could not determine amino acid order to calculate dependencies first.'
				' Make sure there are no cyclical pathways for amino acids that can degrade.')

		for amino_acid in ordered_aa_ids:
			data = self.aa_synthesis_pathways[amino_acid]
			enzymes = data['enzymes']
			enzyme_counts = basal_container.counts(enzymes).sum()

			aa_conc = minimal_conc[amino_acid]
			if data['ki'] is None:
				ki = np.inf * units.mol / units.L
			else:
				# Get largest dynamic range possible given the range of measured KIs
				lower_limit, upper_limit = data['ki']
				if aa_conc < lower_limit:
					ki = lower_limit
				elif aa_conc > upper_limit:
					ki = upper_limit
				else:
					ki = aa_conc
			upstream_aa = [aa for aa in data['upstream']]
			km_conc = METABOLITE_CONCENTRATION_UNITS * np.array([
				minimal_conc[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
				for aa in upstream_aa
				])
			kms_upstream = data['km, upstream']
			kms = METABOLITE_CONCENTRATION_UNITS * np.array([
				kms_upstream.get(aa, minimal_conc[aa]).asNumber(METABOLITE_CONCENTRATION_UNITS)
				for aa in upstream_aa
				])  # TODO: better way to fill in this missing data
			if data['reverse']:
				if np.isnan(data['km, reverse'].asNumber()):
					km_reverse = minimal_conc[amino_acid] * 10  # TODO: better way to fill in this missing data
				else:
					km_reverse = data['km, reverse']
			else:
				km_reverse = np.inf * units.mol / units.L
			km_degradation = data['km, degradation']

			total_supply = supply[amino_acid]
			for aa, stoich in data['downstream'].items():
				total_supply += stoich * supply[aa]
				total_supply += stoich * degradation_rates.get(aa, 0 / units.min)

			# Make required adjustments in order to get positive kcats and import rates
			for parameter, factor in self.aa_synthesis_pathway_adjustments.get(amino_acid, {}).items():
				if parameter == 'ki':
					ki *= factor
				elif parameter == 'km_degradation':
					km_degradation *= factor
				elif parameter == 'km_reverse':
					km_reverse *= factor
				elif parameter == 'kms':
					kms *= factor
				else:
					raise ValueError(f'Unexpected parameter adjustment ({parameter}) for {amino_acid}.')

			# Calculate kcat value to ensure sufficient supply to double
			kcat = total_supply / (enzyme_counts * (1 / (1 + aa_conc / ki) * np.prod(1 / (1 + kms / km_conc)) - 1 / (1 + km_reverse / aa_conc) - 1 / (1 + km_degradation / aa_conc)))
			data['kcat'] = kcat

			aa_enzymes += enzymes
			enzyme_to_aa += [amino_acid] * len(enzymes)
			aa_kcats[amino_acid] = kcat.asNumber(K_CAT_UNITS)
			aa_kis[amino_acid] = ki.asNumber(METABOLITE_CONCENTRATION_UNITS)
			upstream_aas_for_km[amino_acid] = upstream_aa
			upstream_aas[amino_acid] = data['upstream']
			reverse_aas[amino_acid] = data['reverse']
			aa_upstream_kms[amino_acid] = kms.asNumber(METABOLITE_CONCENTRATION_UNITS)
			aa_reverse_kms[amino_acid] = km_reverse.asNumber(METABOLITE_CONCENTRATION_UNITS)
			aa_degradation_kms[amino_acid] = km_degradation.asNumber(METABOLITE_CONCENTRATION_UNITS)
			degradation_rates[amino_acid] = kcat * enzyme_counts / (1 + km_degradation / aa_conc)

		self.aa_enzymes = np.unique(aa_enzymes)
		self.aa_kcats = np.array([aa_kcats[aa] for aa in aa_ids])
		self.aa_kis = np.array([aa_kis[aa] for aa in aa_ids])
		self.aa_upstream_kms = [aa_upstream_kms[aa] for aa in aa_ids]
		self.aa_reverse_kms = np.array([aa_reverse_kms[aa] for aa in aa_ids])
		self.aa_degradation_kms = np.array([aa_degradation_kms[aa] for aa in aa_ids])

		# Convert aa_conc to array with upstream aa_conc via indexing (aa_conc[self.aa_upstream_mapping])
		aa_to_index = {aa: i for i, aa in enumerate(aa_ids)}

		# TODO: better way of handling this that is efficient computationally
		self.aa_to_index = aa_to_index
		self.aa_upstream_aas = [upstream_aas_for_km[aa] for aa in aa_ids]

		# Convert enzyme counts to an amino acid basis via dot product (counts @ self.enzyme_to_amino_acid)
		self.enzyme_to_amino_acid = np.zeros((len(self.aa_enzymes), len(aa_ids)))
		enzyme_mapping = {e: i for i, e in enumerate(self.aa_enzymes)}
		aa_mapping = {a: i for i, a in enumerate(aa_ids)}
		for enzyme, aa in zip(aa_enzymes, enzyme_to_aa):
			self.enzyme_to_amino_acid[enzyme_mapping[enzyme], aa_mapping[aa]] = 1

		# Convert individual supply calculations to overall supply based on dependencies
		# via dot product (self.aa_forward_stoich @ supply)
		# TODO: check for loops (eg ser dependent on glt, glt dependent on ser)
		# TODO: check this makes sense with new format
		self.aa_forward_stoich = np.eye(len(aa_ids))
		for aa, upstream in upstream_aas.items():
			for upstream_aa, stoich in upstream.items():
				self.aa_forward_stoich[aa_to_index[upstream_aa], aa_to_index[aa]] = -stoich
		self.aa_reverse_stoich = np.eye(len(aa_ids))
		for aa, reverse in reverse_aas.items():
			for reverse_aa, stoich in reverse.items():
				self.aa_reverse_stoich[aa_to_index[reverse_aa], aa_to_index[aa]] = -stoich

		# Calculate import rates to match supply in amino acid conditions
		with_aa_rates = (
			sim_data.translation_supply_rate['minimal_plus_amino_acids']
			* cell_specs['with_aa']['avgCellDryMassInit'] * sim_data.constants.n_avogadro
			).asNumber(K_CAT_UNITS)
		with_aa_supply = {
			aa: rate
			for aa, rate in zip(sim_data.molecule_groups.amino_acids, with_aa_rates)
			}
		enzyme_counts = with_aa_container.counts(self.aa_enzymes)
		aa_conc = units.mol / units.L * np.array([
			conc('minimal_plus_amino_acids')[aa].asNumber(units.mol/units.L)
			for aa in aa_ids
			])

		supply = np.array([with_aa_supply[aa] for aa in aa_ids])  # TODO: check that this is ok and do not need other adjustments for downstream
		synthesis, _, _ = self.amino_acid_synthesis(enzyme_counts, aa_conc)
		self.specific_import_rates = (supply - synthesis) / cell_specs['with_aa']['avgCellDryMassInit'].asNumber(DRY_MASS_UNITS)
		self.specific_import_rates[aa_ids.index('CYS[c]')] = 0  # Does not have direct import

		# Concentrations for reference in analysis plot
		conversion = sim_data.constants.cell_density / sim_data.constants.n_avogadro * sim_data.mass.cell_dry_mass_fraction
		basal_counts = basal_container.counts(self.aa_enzymes)
		self.aa_supply_enzyme_conc_with_aa = conversion * enzyme_counts / cell_specs['with_aa']['avgCellDryMassInit']
		self.aa_supply_enzyme_conc_basal = conversion * basal_counts / cell_specs['basal']['avgCellDryMassInit']

		# Check calculations that could end up negative
		neg_idx = np.where(self.aa_kcats < 0)[0]
		if len(neg_idx):
			aas = ', '.join([aa_ids[idx] for idx in neg_idx])
			print(f'{self.aa_kcats = }')
			raise ValueError(f'kcat value was determined to be negative for {aas}.'
				' Check input parameters like KM and KI or the concentration.')
		neg_idx = np.where(self.specific_import_rates < 0)[0]
		if len(neg_idx):
			aas = ', '.join([aa_ids[idx] for idx in neg_idx])
			print(f'{self.specific_import_rates = }')
			raise ValueError(f'Import rate was determined to be negative for {aas}.'
				' Check input parameters like supply and synthesis or enzyme expression.')

	def amino_acid_synthesis(self, enzyme_counts: np.ndarray, aa_conc: units.Unum):
		"""
		Calculate the net rate of synthesis for amino acid pathways (can be
		negative with reverse reactions).

		Args:
			enzyme_counts: counts for each enzyme accounted for in pathways
			aa_conc: concentrations of each amino acid with mol/volume units

		Returns:
			synthesis: net rate of synthesis for each amino acid pathway.
				array is unitless but represents counts of amino acid per second
		"""

		# Convert to appropraite arrays
		aa_conc = aa_conc.asNumber(METABOLITE_CONCENTRATION_UNITS)
		counts_per_aa = enzyme_counts @ self.enzyme_to_amino_acid

		# TODO: more efficient way of doing this
		km_saturation = np.array([np.product([1 / (1 + km / aa_conc[self.aa_to_index[aa]]) for km, aa in zip(kms, aas)]) for kms, aas in zip(self.aa_upstream_kms, self.aa_upstream_aas)])

		# Determine saturation fraction for reactions
		forward_fraction = 1 / (1 + aa_conc / self.aa_kis) * km_saturation
		reverse_fraction = 1 / (1 + self.aa_reverse_kms / aa_conc)
		loss_fraction = 1 / (1 + self.aa_degradation_kms / aa_conc)
		fraction = forward_fraction - reverse_fraction - loss_fraction

		# Calculate synthesis rate
		synthesis = (
			self.aa_forward_stoich @ (self.aa_kcats * counts_per_aa * forward_fraction)
			- self.aa_reverse_stoich @ (self.aa_kcats * counts_per_aa * reverse_fraction)
			- self.aa_kcats * counts_per_aa * loss_fraction
		)

		return synthesis, counts_per_aa, fraction

	def amino_acid_export(self, aa_transporters_counts: np.ndarray, aa_conc: units.Unum, mechanistic_uptake: bool):
		"""
		Calculate the rate of amino acid export.

		Args:
			aa_transporters_counts: counts of each transporter
			aa_conc: concentrations of each amino acid with mol/volume units
			mechanisitc_uptake: if true, the uptake is calculated based on transporters

		Returns:
			rate of export for each amino acid. array is unitless but
				represents counts of amino acid per second
		"""

		if mechanistic_uptake:
			# Export based on mechanistic model
			aa_conc = aa_conc.asNumber(METABOLITE_CONCENTRATION_UNITS)
			trans_counts_per_aa = self.aa_to_export_transporters_matrix @ aa_transporters_counts
			export_rates = self.export_kcats_per_aa * trans_counts_per_aa / (1 + self.aa_export_kms / aa_conc)
		else:
			# Export is lumped with specific uptake rates in amino_acid_import
			# and not dependent on internal amino acid concentrations or
			# explicitly considered here
			export_rates = np.zeros(len(aa_conc))

		return export_rates

	def amino_acid_import(self, aa_in_media: np.ndarray, dry_mass: units.Unum, aa_transporters_counts: np.ndarray, mechanisitc_uptake: bool):
		"""
		Calculate the rate of amino acid uptake.

		Args:
			aa_in_media: bool for each amino acid being present in current media
			dry_mass: current dry mass of the cell, with mass units
			aa_transporters_counts: counts of each transporter
			mechanisitc_uptake: if true, the uptake is calculated based on transporters

		Returns:
			rate of uptake for each amino acid. array is unitless but
				represents counts of amino acid per second
		"""

		if mechanisitc_uptake:
			# Uptake based on mechanistic model
			counts_per_aa = self.aa_to_transporters_matrix @ aa_transporters_counts
			import_rates = self.uptake_kcats_per_aa * counts_per_aa
		else:
			import_rates = self.specific_import_rates * dry_mass.asNumber(DRY_MASS_UNITS)

		return import_rates * aa_in_media

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
		compartment_ids_to_abbreviations = {
			comp['id']: comp['abbrev'] for comp in raw_data.compartments
			}
		compartment_ids_to_abbreviations.update(
			UNDEFINED_COMPARTMENT_IDS_TO_ABBREVS)

		valid_directions = {'L2R', 'R2L', 'BOTH'}
		forward_directions = {'L2R', 'BOTH'}
		reverse_directions = {'R2L', 'BOTH'}

		metabolite_ids = {met['id'] for met in cast(Any, raw_data).metabolites}

		# Build mapping from each complexation subunit to all downstream
		# complexes containing the subunit, including itself
		# Start by building mappings from subunits to complexes that are
		# directly formed from the subunit through a single reaction
		subunit_id_to_parent_complexes = {} # type: Dict[str, List[str]]

		for comp_reaction in itertools.chain(
				cast(Any, raw_data).complexation_reactions,
				cast(Any, raw_data).equilibrium_reactions):
			complex_id = None

			# Find ID of complex
			for mol_id, coeff in comp_reaction['stoichiometry'].items():
				if coeff > 0:
					complex_id = mol_id
					break

			assert complex_id is not None

			# Map each subunit to found complex
			for mol_id, coeff in comp_reaction['stoichiometry'].items():
				if mol_id == complex_id or mol_id in metabolite_ids:
					continue
				elif mol_id in subunit_id_to_parent_complexes:
					subunit_id_to_parent_complexes[mol_id].append(complex_id)
				else:
					subunit_id_to_parent_complexes[mol_id] = [complex_id]

		# Recursive function that returns a list of all downstream complexes
		# containing the given subunit, including itself
		def get_all_complexes(subunit_id):
			all_downstream_complex_ids = [subunit_id]

			if subunit_id not in subunit_id_to_parent_complexes:
				return all_downstream_complex_ids

			# Get downstream complexes of all parent complexes
			for parent_complex_id in subunit_id_to_parent_complexes[subunit_id]:
				all_downstream_complex_ids.extend(get_all_complexes(parent_complex_id))

			# Remove duplicates
			return sorted(set(all_downstream_complex_ids))

		subunit_id_to_all_downstream_complexes = {
			subunit_id: get_all_complexes(subunit_id)
			for subunit_id in subunit_id_to_parent_complexes.keys()
			}

		# Initialize variables to store reaction information
		reaction_stoich = {}
		reversible_reactions = []
		reaction_catalysts = {}

		# Load and parse reaction information from raw_data
		for reaction in cast(Any, raw_data).metabolic_reactions:
			reaction_id = reaction["id"]
			stoich = reaction["stoichiometry"]
			direction = reaction["direction"]

			if len(stoich) <= 1:
				raise Exception("Invalid biochemical reaction: {}, {}".format(reaction_id, stoich))

			if direction not in valid_directions:
				raise InvalidReactionDirectionError(
					f'The direction {direction} given for reaction {reaction_id} is invalid.')

			forward = direction in forward_directions
			reverse = direction in reverse_directions

			def convert_compartment_tags(met_id):
				new_met_id = met_id

				for comp_id, comp_abbrev in compartment_ids_to_abbreviations.items():
					new_met_id = new_met_id.replace(
						f'[{comp_id}]', f'[{comp_abbrev}]')

				return new_met_id

			# All protein complexes that contain an enzyme subunit are assumed
			# to retain the enzyme's catalytic activity
			catalysts_for_this_rxn = []
			all_potential_catalysts = []
			for catalyst in reaction["catalyzed_by"]:
				all_potential_catalysts.extend(
					subunit_id_to_all_downstream_complexes.get(catalyst, [catalyst]))

			for catalyst in sorted(set(all_potential_catalysts)):
				if sim_data.getter.is_valid_molecule(catalyst):
					catalysts_with_loc = catalyst + sim_data.getter.get_compartment_tag(catalyst)
					catalysts_for_this_rxn.append(catalysts_with_loc)
				# If we don't have the catalyst in our reconstruction, drop it
				else:
					if VERBOSE:
						print(
							'Skipping catalyst {} for {} since it is not in the model'
							.format(catalyst, reaction_id)
							)

			if forward:
				reaction_stoich[reaction_id] = {
					convert_compartment_tags(moleculeID): stoichCoeff
					for moleculeID, stoichCoeff in stoich.items()
					}
				if len(catalysts_for_this_rxn) > 0:
					reaction_catalysts[reaction_id] = catalysts_for_this_rxn

			if reverse:
				reverse_reaction_id = REVERSE_REACTION_ID.format(reaction_id)
				reaction_stoich[reverse_reaction_id] = {
					convert_compartment_tags(moleculeID): -stoichCoeff
					for moleculeID, stoichCoeff in stoich.items()
					}
				if len(catalysts_for_this_rxn) > 0:
					reaction_catalysts[reverse_reaction_id] = list(catalysts_for_this_rxn)

			if forward and reverse:
				reversible_reactions.append(reaction_id)

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

	def _is_transport_rxn(self, stoich):
		# type: (Dict[str, int]) -> bool
		"""
		Determines if the metabolic reaction with a given stoichiometry is a
		transport reactions that transports metabolites between different
		compartments. A metabolic reaction is considered to be a transport
		reaction if the substrate set and the product share the same metabolite
		tagged into different compartments.

		Args:
			stoich: Stoichiometry of the metabolic reaction
				{
				metabolite ID (str): stoichiometric coefficient (int)
				}

		Returns:
			is_transport_rxn (bool): True if the reaction with the given
			stoichiometry is a transport reaction
		"""
		is_transport_rxn = False

		# Get IDs of all substrates and products
		substrates, products = [], []
		for mol_id, coeff in stoich.items():
			if coeff < 0:
				substrates.append(mol_id)
			else:
				products.append(mol_id)

		# Get mapping from IDs to IDs without compartments
		substrates_tagged_to_no_tag = {
			mol_id: re.sub(r"\[.*]", "", mol_id) for mol_id in substrates}
		products_tagged_to_no_tag = {
			mol_id: re.sub(r"\[.*]", "", mol_id) for mol_id in products}

		overlap_no_tag = (
			set(substrates_tagged_to_no_tag.values())
			& set(products_tagged_to_no_tag.values()))

		for mol_id_no_tag in list(overlap_no_tag):
			substrates_tagged = [
				mol_tagged for mol_tagged in substrates
				if substrates_tagged_to_no_tag[mol_tagged] == mol_id_no_tag]
			products_tagged = [
				mol_tagged for mol_tagged in products
				if products_tagged_to_no_tag[mol_tagged] == mol_id_no_tag]

			overlap_tagged = set(substrates_tagged) & set(products_tagged)

			# Tag reaction as a transport reaction if there is no overlap
			# between those substrates and products with locations included
			if len(overlap_tagged) == 0:
				is_transport_rxn = True
				break

		return is_transport_rxn


# Class used to update metabolite concentrations based on the current nutrient conditions
class ConcentrationUpdates(object):
	def __init__(self, concDict, relative_changes, equilibriumReactions, exchange_data_dict, all_metabolite_ids, linked_metabolites):
		self.units = units.getUnit(list(concDict.values())[0])
		self.default_concentrations_dict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)
		self.exchange_fluxes = self._exchange_flux_present(exchange_data_dict)
		self.relative_changes = relative_changes
		self._all_metabolite_ids = all_metabolite_ids
		self.linked_metabolites = linked_metabolites

		# factor of internal amino acid increase if amino acids present in nutrients
		self.molecule_scale_factors = {
			"L-ALPHA-ALANINE[c]": 2.,
			"ARG[c]": 2.,
			"ASN[c]": 2.,
			"L-ASPARTATE[c]": 2.,
			"CYS[c]": 2.,
			"GLT[c]": 2.,
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

		self.molecule_set_amounts = self._add_molecule_amounts(equilibriumReactions, self.default_concentrations_dict)

	# return adjustments to concDict based on nutrient conditions
	def concentrations_based_on_nutrients(self, media_id=None, imports=None, conversion_units=None):
		# type: (Optional[str], Optional[Set[str]], Optional[units.Unum]) -> Dict[str, Any]
		if conversion_units:
			conversion = self.units.asNumber(conversion_units)
		else:
			conversion = self.units

		if imports is None and media_id is not None:
			imports = self.exchange_fluxes[media_id]

		concentrationsDict = self.default_concentrations_dict.copy()

		metaboliteTargetIds = sorted(concentrationsDict.keys())
		concentrations = conversion * np.array([concentrationsDict[k] for k in metaboliteTargetIds])
		concDict = dict(zip(metaboliteTargetIds, concentrations))

		if imports is not None:
			# For faster conversions than .asNumber(conversion_units) for each setAmount
			if conversion_units:
				conversion_to_no_units = conversion_units.asUnit(self.units)

			# Adjust for measured concentration changes in different media
			if media_id in self.relative_changes:
				for mol_id, conc_change in self.relative_changes[media_id].items():
					if mol_id in concDict:
						concDict[mol_id]  *= conc_change

			for moleculeName, setAmount in self.molecule_set_amounts.items():
				if ((moleculeName in imports and (moleculeName[:-3] + "[c]" not in self.molecule_scale_factors or moleculeName == "L-SELENOCYSTEINE[c]"))
						or (moleculeName in self.molecule_scale_factors and moleculeName[:-3] + "[p]" in imports)):
					if conversion_units:
						setAmount = (setAmount / conversion_to_no_units).asNumber()
					concDict[moleculeName] = setAmount

		for met, linked in self.linked_metabolites.items():
			concDict[met] = concDict[linked['lead']] * linked['ratio']

		return concDict

	def _exchange_flux_present(self, exchange_data):
		# type: (Dict[str, Any]) -> Dict[str, Set[str]]
		"""
		Caches the presence of exchanges in each media condition based on
		exchange_data to set concentrations in concentrations_based_on_nutrients().

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

	def _add_molecule_amounts(self, equilibriumReactions, concDict):
		moleculeSetAmounts = {}
		for reaction in equilibriumReactions:
			# We only want to do this for species with standard Michaelis-Menten kinetics initially
			if len(reaction["stoichiometry"]) != 3:
				continue

			moleculeName = [
				mol_id for mol_id in reaction["stoichiometry"].keys()
				if mol_id in self._all_metabolite_ids][0]
			amountToSet = 1e-4
			moleculeSetAmounts[moleculeName + "[p]"] = amountToSet * self.units
			moleculeSetAmounts[moleculeName + "[c]"] = amountToSet * self.units

		for moleculeName, scaleFactor in six.viewitems(self.molecule_scale_factors):
			moleculeSetAmounts[moleculeName] = scaleFactor * concDict[moleculeName] * self.units
		return moleculeSetAmounts
