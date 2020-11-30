"""
SimulationData getter functions
"""

from __future__ import absolute_import, division, print_function

import itertools
import re
from typing import Any, List, Union

from Bio.Seq import Seq
import numpy as np

from reconstruction.ecoli.dataclasses.molecule_groups import POLYMERIZED_FRAGMENT_PREFIX
from wholecell.utils import units


class UnknownMolecularWeightError(Exception):
	pass

class UnknownModifiedRnaError(Exception):
	pass

class InvalidProteinComplexError(Exception):
	pass


class GetterFunctions(object):
	""" getterFunctions """

	def __init__(self, raw_data, sim_data):
		self._n_submass_indexes = len(sim_data.submass_name_to_index)
		self._submass_name_to_index = sim_data.submass_name_to_index

		self._build_sequences(raw_data)
		self._build_all_masses(raw_data, sim_data)
		self._build_compartments(raw_data, sim_data)

	def get_sequences(self, ids):
		# type: (Union[List[str], np.ndarray]) -> List[str]
		"""
		Return a list of sequences of the molecules with the given IDs.
		"""
		assert isinstance(ids, (list, np.ndarray))
		return [self._sequences[mol_id] for mol_id in ids]

	def get_mass(self, mol_id):
		# type: (str) -> Any
		"""
		Return the total mass of the molecule with a given ID.
		"""
		assert isinstance(mol_id, str)
		return self._mass_units * self._all_total_masses[self._compartment_tag.sub('', mol_id)]

	def get_masses(self, mol_ids):
		# type: (Union[List, np.ndarray]) -> Any
		"""
		Return an array of total masses of the molecules with the given IDs.
		"""
		assert isinstance(mol_ids, (list, np.ndarray))
		masses = [
			self._all_total_masses[self._compartment_tag.sub('', mol_id)]
			for mol_id in mol_ids]
		return self._mass_units * np.array(masses)

	def get_submass_array(self, mol_id):
		# type: (str) -> Any
		"""
		Return the submass array of the molecule with a given ID.
		"""
		assert isinstance(mol_id, str)
		return self._mass_units * self._all_submass_arrays[self._compartment_tag.sub('', mol_id)]

	def get_compartment(self, mol_id):
		# type: (str) -> List[List[str]]
		"""
		Returns the list of one-letter codes for the compartments that the
		molecule with the given ID can exist in.
		"""
		assert isinstance(mol_id, str)
		return self._all_compartments[mol_id]

	def get_compartments(self, ids):
		# type: (Union[List[str], np.ndarray]) -> List[List[str]]
		"""
		Returns a list of the list of one-letter codes for the compartments
		that each of the molecules with the given IDs can exist in.
		"""
		assert isinstance(ids, (list, np.ndarray))
		return [self._all_compartments[x] for x in ids]

	def get_compartment_tag(self, id_):
		# type: (str) -> str
		"""
		Look up a molecule id and return a compartment suffix tag like '[c]'.
		"""
		return f'[{self._all_compartments[id_][0]}]'

	def check_valid_molecule(self, mol_id):
		# type: (str) -> bool
		"""
		Returns True if the molecule with the given ID is a valid molecule (has
		both a submass array and a compartment tag).
		"""
		return mol_id in self._all_submass_arrays and mol_id in self._all_compartments

	def get_all_valid_molecules(self):
		# type: () -> list
		"""
		Returns a list of all molecule IDs with assigned submass arrays and
		compartments.
		"""
		return sorted(self._all_submass_arrays.keys() & self._all_compartments.keys())

	def _build_sequences(self, raw_data):
		"""
		Builds sequences of RNAs and proteins.
		"""
		self._sequences = {}
		self._build_rna_sequences(raw_data)
		self._build_protein_sequences(raw_data)

	def _build_rna_sequences(self, raw_data):
		"""
		Builds nucleotide sequences of each RNA using the genome sequence and
		the transcription start sites and lengths of the corresponding gene.
		"""
		# Get index of gene corresponding to each RNA
		rna_id_to_gene_index = {gene['rna_id']: i
			for i, gene in enumerate(raw_data.genes)}

		# Get RNA lengths from gene data
		gene_lengths = [gene['length'] for gene in raw_data.genes]

		# Get list of coordinates and directions for each gene
		coordinate_list = [gene["coordinate"] for gene in raw_data.genes]
		direction_list = [gene["direction"] for gene in raw_data.genes]

		# Get RNA sequence from genome sequence
		genome_sequence = raw_data.genome_sequence

		for rna in raw_data.rnas:
			gene_index = rna_id_to_gene_index[rna['id']]
			coordinate = coordinate_list[gene_index]

			# Parse genome sequence to get RNA sequence
			if direction_list[gene_index] == '+':
				seq = genome_sequence[coordinate:coordinate + gene_lengths[gene_index]].transcribe()
			else:
				seq = genome_sequence[coordinate - gene_lengths[gene_index] + 1:coordinate + 1].reverse_complement().transcribe()

			self._sequences[rna['id']] = seq

	def _build_protein_sequences(self, raw_data):
		"""
		Builds the amino acid sequences of each protein monomer using sequences
		in raw_data.
		"""
		for protein in raw_data.proteins:
			self._sequences[protein['id']] = Seq(protein['seq'])

	def _build_submass_array(self, mw, submass_name):
		# type: (float, str) -> np.ndarray
		"""
		Converts a scalar molecular weight value to an array of submasses,
		given the name of the submass category that the molecule belongs to.
		"""
		mw_array = np.zeros(self._n_submass_indexes)
		mw_array[self._submass_name_to_index[submass_name]] = mw
		return mw_array

	def _build_all_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights keyed with the IDs of molecules.
		Depending on the type of the molecule, weights are either pulled
		directly from raw data, or dynamically calculated from the existing
		data.
		"""
		self._all_submass_arrays = {}

		self._all_submass_arrays.update({
			met['id']: (self._build_submass_array(met['mw'], 'metabolite')
			if met['id'] != sim_data.molecule_ids.water[:-3]
			else self._build_submass_array(met['mw'], 'water'))
			for met in raw_data.metabolites
			})

		# These updates can be dependent on metabolite masses
		self._all_submass_arrays.update(self._build_polymerized_subunit_masses(sim_data))
		self._all_submass_arrays.update(self._build_rna_masses(raw_data, sim_data))
		self._all_submass_arrays.update(self._build_protein_masses(raw_data, sim_data))
		self._all_submass_arrays.update(self._build_full_chromosome_mass(raw_data, sim_data))

		# These updates can be dependent on the masses calculated above
		self._all_submass_arrays.update(self._build_modified_rna_masses(raw_data))
		self._all_submass_arrays.update(self._build_protein_complex_masses(raw_data))
		self._all_submass_arrays.update(
			{x['id']: np.array(x['mw']) for x in raw_data.modified_proteins}
			)

		self._mass_units = units.g / units.mol
		self._compartment_tag = re.compile(r'\[[a-z]\]')

		# Build dictionary of total masses of each molecule
		self._all_total_masses = {
			mol_id: mass_array.sum()
			for (mol_id, mass_array) in self._all_submass_arrays.items()
			}

	def _build_polymerized_subunit_masses(self, sim_data):
		"""
		Builds dictionary of molecular weights of polymerized subunits (NTPs,
		dNTPs, and amino acids) by subtracting the end weights from the weights
		of original subunits.
		"""
		polymerized_subunit_masses = {}

		def add_subunit_masses(subunit_ids, end_group_id):
			"""
			Add subunit masses to dictionary given the IDs of the original
			subunits and the end group.
			"""
			# Subtract mw of end group from each polymer subunit
			end_group_mw = self._all_submass_arrays[end_group_id[:-3]].sum()
			polymerized_subunit_mws = [
				self._all_submass_arrays[met_id[:-3]].sum() - end_group_mw
				for met_id in subunit_ids]

			# Add to dictionary with prefixed ID
			polymerized_subunit_masses.update({
				POLYMERIZED_FRAGMENT_PREFIX + subunit_id[:-3]: self._build_submass_array(mw, 'metabolite')
					for subunit_id, mw in zip(subunit_ids, polymerized_subunit_mws)
				})

		add_subunit_masses(
			sim_data.molecule_groups.ntps, sim_data.molecule_ids.ppi)
		add_subunit_masses(
			sim_data.molecule_groups.dntps, sim_data.molecule_ids.ppi)
		add_subunit_masses(
			sim_data.molecule_groups.amino_acids, sim_data.molecule_ids.water)

		return polymerized_subunit_masses

	def _build_rna_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights of RNAs keyed with the RNA IDs.
		Molecular weights are calculated from the RNA sequence and the weights
		of polymerized NTPs.
		"""
		# Get RNA nucleotide compositions
		rna_seqs = self.get_sequences([rna['id'] for rna in raw_data.rnas])
		nt_counts = []
		for seq in rna_seqs:
			nt_counts.append(
				[seq.count(letter) for letter in sim_data.ntp_code_to_id_ordered.keys()])
		nt_counts = np.array(nt_counts)

		# Calculate molecular weights
		ppi_mw = self._all_submass_arrays[sim_data.molecule_ids.ppi[:-3]].sum()
		polymerized_ntp_mws = np.array([
			self._all_submass_arrays[met_id[:-3]].sum() for met_id in sim_data.molecule_groups.polymerized_ntps
			])

		mws = nt_counts.dot(polymerized_ntp_mws) + ppi_mw  # Add end weight

		return {rna['id']: self._build_submass_array(mw, rna['type']) for (rna, mw) in zip(raw_data.rnas, mws)}

	def _build_protein_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights of protein monomers keyed with
		the protein IDs. Molecular weights are calculated from the protein
		sequence and the weights of polymerized amino acids.
		"""
		# Get protein amino acid compositions
		protein_seqs = self.get_sequences(
			[protein['id'] for protein in raw_data.proteins])
		aa_counts = []
		for seq in protein_seqs:
			aa_counts.append(
				[seq.count(letter) for letter in sim_data.amino_acid_code_to_id_ordered.keys()])
		aa_counts = np.array(aa_counts)

		# Calculate molecular weights
		water_mw = self._all_submass_arrays[sim_data.molecule_ids.water[:-3]].sum()
		polymerized_aa_mws = np.array(
			[self._all_submass_arrays[met_id[:-3]].sum() for met_id in sim_data.molecule_groups.polymerized_amino_acids]
			)
		mws = aa_counts.dot(polymerized_aa_mws) + water_mw  # Add end weight

		return {protein['id']: self._build_submass_array(mw, 'protein') for (protein, mw) in zip(raw_data.proteins, mws)}

	def _build_full_chromosome_mass(self, raw_data, sim_data):
		"""
		Calculates the mass of the full chromosome from its sequence and the
		weigths of polymerized dNTPs.
		"""
		# Get chromosome dNTP compositions
		chromosome_seq = raw_data.genome_sequence
		forward_strand_nt_counts = np.array([
			chromosome_seq.count(letter)
			for letter in sim_data.dntp_code_to_id_ordered.keys()
			])
		reverse_strand_nt_counts = np.array([
			chromosome_seq.reverse_complement().count(letter)
			for letter in sim_data.dntp_code_to_id_ordered.keys()
			])

		# Calculate molecular weight
		polymerized_dntp_mws = np.array([
			self._all_submass_arrays[met_id[:-3]].sum() for met_id in sim_data.molecule_groups.polymerized_dntps
			])
		mw = float(np.dot(
			forward_strand_nt_counts + reverse_strand_nt_counts,
			polymerized_dntp_mws))

		return {sim_data.molecule_ids.full_chromosome[:-3]: self._build_submass_array(mw, 'DNA')}

	def _build_modified_rna_masses(self, raw_data):
		"""
		Builds dictionary of molecular weights of modified RNAs keyed with
		the molecule IDs. Molecular weights are calculated from the
		stoichiometries of the modification reactions.
		"""
		modified_rna_masses = {}

		# Get all modified form IDs from RNA data
		all_modified_rna_ids = {
			modified_rna_id for rna in raw_data.rnas
			for modified_rna_id in rna['modified_forms']}

		# Get IDs of modification reactions that should be removed
		removed_rna_modification_reaction_ids = {
			rxn['id'] for rxn in raw_data.rna_modification_reactions_removed
			}

		# Loop through each modification reaction
		for rxn in raw_data.rna_modification_reactions:
			if rxn['id'] in removed_rna_modification_reaction_ids:
				continue

			# Find molecule IDs whose masses are unknown
			unknown_mol_ids = [
				x['molecule'] for x in rxn['stoichiometry']
				if x['molecule'] not in self._all_submass_arrays]
			# The only molecule whose mass is unknown should be the modified
			# RNA molecule
			if len(unknown_mol_ids) != 1:
				raise UnknownMolecularWeightError(
					'RNA modification reaction %s has two or more reactants or products with unknown molecular weights: %s' % (rxn['id'], unknown_mol_ids))

			# All modified RNAs must be assigned to a specific RNA
			modified_rna_id = unknown_mol_ids[0]
			if modified_rna_id not in all_modified_rna_ids:
				raise UnknownModifiedRnaError(
					'The modified RNA molecule %s was not assigned to any RNAs.' % (modified_rna_id, ))

			# Calculate mw of the modified RNA from the stoichiometry
			mw_sum = np.zeros(self._n_submass_indexes)
			for entry in rxn['stoichiometry']:
				if entry['molecule'] == modified_rna_id:
					modified_rna_coeff = entry['coeff']
				else:
					mw_sum += entry['coeff']*self._all_submass_arrays[entry['molecule']]

			modified_rna_masses[modified_rna_id] = -mw_sum/modified_rna_coeff

		return modified_rna_masses

	def _build_protein_complex_masses(self, raw_data):
		"""
		Builds dictionary of molecular weights of protein complexes keyed with
		the molecule IDs. Molecular weights are calculated from the
		stoichiometries of the complexation/equilibrium reactions. For
		complexes whose subunits are also complexes, the molecular weights are
		calculated recursively.
		"""
		protein_complex_masses = {}

		# Get IDs of complexation/equilibrium reactions that should be removed
		removed_reaction_ids = {
			rxn['id'] for rxn in itertools.chain(
				raw_data.complexation_reactions_removed,
				raw_data.equilibrium_reactions_removed)
			}

		# Build mapping from complex ID to its subunit stoichiometry
		complex_id_to_stoich = {}

		for rxn in itertools.chain(
				raw_data.complexation_reactions, raw_data.equilibrium_reactions):
			# Skip removed reactions
			if rxn['id'] in removed_reaction_ids:
				continue

			# Get the ID of the complex and the stoichiometry of subunits
			complex_ids = []
			subunit_stoich = {}

			for mol in rxn['stoichiometry']:
				if mol['coeff'] == 1:
					complex_ids.append(mol['molecule'])
				elif mol['coeff'] < 0:
					subunit_stoich.update({mol['molecule']: -mol['coeff']})

				# Each molecule should either be a complex with the coefficient
				# 1 or a subunit with a negative coefficient
				else:
					raise InvalidProteinComplexError(
						'Reaction %s contains an invalid molecule %s with stoichiometric coefficient %.1f.'
						% (rxn['id'], mol['molecule'], mol['coeff']))

			# Each reaction should only produce a single protein complex
			if len(complex_ids) != 1:
				raise InvalidProteinComplexError(
					'Reaction %s does not produce a single protein complex.'
					% (rxn['id'], ))

			complex_id_to_stoich[complex_ids[0]] = subunit_stoich

		def get_mw(mol_id):
			"""
			Recursive function that returns the molecular weight given a
			molecule ID.
			"""
			# Look for molecule ID in existing dictionaries, and return
			# molecular weight array if found
			if mol_id in self._all_submass_arrays:
				return self._all_submass_arrays[mol_id]
			elif mol_id in protein_complex_masses:
				return protein_complex_masses[mol_id]

			# If molecule ID does not exist, recursively add up the molecular
			# weights of the molecule's subunits
			else:
				# Each complex molecule should have a corresponding subunit
				# stoichiometry
				if mol_id not in complex_id_to_stoich:
					raise InvalidProteinComplexError(
						'Complex %s is not being produced by any copmlexation or equilibrium reaction.'
						% (mol_id, )
						)
				mw_sum = np.zeros(self._n_submass_indexes)

				for subunit_id, coeff in complex_id_to_stoich[mol_id].items():
					mw_sum += coeff * get_mw(subunit_id)

				return mw_sum

		# Get molecular weights of each protein complex
		for complex_id in complex_id_to_stoich.keys():
			protein_complex_masses.update({complex_id: get_mw(complex_id)})

		return protein_complex_masses

	def _build_compartments(self, raw_data, sim_data):
		self._all_compartments = {}
		all_compartments = [comp['abbrev'] for comp in raw_data.compartments]

		# RNAs, modified RNAs, and full chromosomes only localize to the
		# cytosol
		self._all_compartments.update({
			rna['id']: ['c'] for rna in raw_data.rnas
			})
		self._all_compartments.update({
			modified_rna_id: ['c'] for rna in raw_data.rnas
			for modified_rna_id in rna['modified_forms']
			if modified_rna_id in self._all_submass_arrays
			})
		self._all_compartments.update({
			sim_data.molecule_ids.full_chromosome[:-3]: ['c']
			})

		# Proteins and modified proteins localize to the single compartment
		# specified in raw data for each protein
		self._all_compartments.update({
			protein['id']: [protein['compartment']] for protein
			in itertools.chain(raw_data.proteins, raw_data.modified_proteins)
			})

		# Metabolites and polymerized subunits can localize to all compartments
		# being modeled
		self._all_compartments.update({
			met['id']: all_compartments for met in raw_data.metabolites
			})
		self._all_compartments.update({
			subunit_id[:-3]: all_compartments
			for subunit_id in sim_data.molecule_groups.polymerized_subunits
			})

		# Sort compartment abbreviations based on complexation priority: when
		# a molecule in a compartment with a higher complexation priority
		# forms a complex with a molecule in a compartment with low priority,
		# the complex is assigned to the higher-priority compartment.
		complexation_priorities = np.array([
			comp['complexation_priority'] for comp in raw_data.compartments])
		all_compartments_sorted = [
			all_compartments[i] for i in np.argsort(complexation_priorities)]

		# Protein complexes are localized based on the compartments of their
		# subunits
		self._all_compartments.update(
			self._build_protein_complex_compartments(raw_data, all_compartments_sorted))


	def _build_protein_complex_compartments(self, raw_data, all_compartments_sorted):
		"""
		Builds dictionary of compartment tags for protein complexes keyed with
		the molecule IDs. Each complex is assigned to the compartment that
		contains a subunit of the complex and has the highest complexation
		priority. Compartment tags of subunits that are also protein complexes
		are determined recursively. Compartments of metabolite subunits are not
		considered since metabolites are assigned to all compartments that are
		being modeled.
		"""
		protein_complex_compartments = {}
		metabolite_ids = {met['id'] for met in raw_data.metabolites}

		# Get IDs of complexation/equilibrium reactions that should be removed
		removed_reaction_ids = {
			rxn['id'] for rxn in itertools.chain(
				raw_data.complexation_reactions_removed,
				raw_data.equilibrium_reactions_removed)
			}

		# Build mapping from complex ID to its subunit IDs
		complex_id_to_subunit_ids = {}

		for rxn in itertools.chain(
				raw_data.complexation_reactions,
				raw_data.equilibrium_reactions):
			# Skip removed reactions
			if rxn['id'] in removed_reaction_ids:
				continue

			# Get the ID of the complex and the stoichiometry of subunits
			complex_id = None
			subunit_ids = []

			for mol in rxn['stoichiometry']:
				if mol['coeff'] == 1:
					complex_id = mol['molecule']
				elif mol['coeff'] < 0:
					subunit_ids.append(mol['molecule'])

			complex_id_to_subunit_ids[complex_id] = subunit_ids

		def get_compartment(mol_id):
			"""
			Recursive function that returns a compartment tag given a molecule ID.
			"""
			# Look for molecule ID in existing dictionaries, and return
			# compartment tag if found
			if mol_id in self._all_compartments:
				return self._all_compartments[mol_id]
			elif mol_id in protein_complex_compartments:
				return protein_complex_compartments[mol_id]

			# If molecule ID does not exist, return the subunit compartment
			# with the highest priority
			else:
				all_subunit_compartments = []
				for subunit_id in complex_id_to_subunit_ids[mol_id]:
					# Skip metabolites (can exist in all compartments)
					if subunit_id in metabolite_ids:
						continue
					else:
						all_subunit_compartments.extend(get_compartment(subunit_id))

				if len(all_subunit_compartments) == 0:
					raise InvalidProteinComplexError(
						'Complex %s is composed of only metabolites.'
						% (mol_id,)
						)

				# Return compartment with highest complexation priority
				for loc in all_compartments_sorted:
					if loc in all_subunit_compartments:
						return [loc]

		# Get compartments of each protein complex
		for complex_id in complex_id_to_subunit_ids.keys():
			protein_complex_compartments.update(
				{complex_id: get_compartment(complex_id)})

		return protein_complex_compartments
