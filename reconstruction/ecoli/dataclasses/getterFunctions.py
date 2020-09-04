"""
SimulationData getter functions

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import absolute_import, division, print_function

import itertools
import re
from typing import List, Union

from Bio.Seq import Seq
import numpy as np

from reconstruction.ecoli.dataclasses.moleculeGroups import POLYMERIZED_FRAGMENT_PREFIX
from wholecell.utils import units

class getterFunctions(object):
	""" getterFunctions """

	def __init__(self, raw_data, sim_data):
		self._build_sequences(raw_data)
		self._build_all_masses(raw_data, sim_data)
		self._build_locations(raw_data, sim_data)

	def get_sequence(self, ids):
		# type: (Union[List[str], np.ndarray]) -> List[str]
		assert isinstance(ids, (list, np.ndarray))
		return [self._sequences[mol_id] for mol_id in ids]

	def getMass(self, ids):
		assert isinstance(ids, (list, np.ndarray))
		masses = [self._all_mass[self._location_tag.sub('', i)] for i in ids]
		return self._mass_units * np.array(masses)

	def getLocation(self, ids):
		# type: (Union[List[str], np.ndarray]) -> List[str]
		assert isinstance(ids, (list, np.ndarray))
		return [self._locationDict[x] for x in ids]

	def get_location_tag(self, id_):
		# type: (str) -> str
		"""Look up a location id and return a location suffix tag like '[c]'."""
		return f'[{self._locationDict[id_][0]}]'

	def check_valid_molecule(self, mol_id):
		return mol_id in self._all_mass and mol_id in self._locationDict

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
		rna_id_to_gene_index = {gene["rnaId"]: i
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

	def _build_all_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights keyed with the IDs of molecules.
		Depending on the type of the molecule, weights are either pulled
		directly from raw data, or dynamically calculated from the existing
		data.
		"""
		self._all_mass = {}

		self._all_mass.update(
			{x['id']: np.sum(x['mw']) for x in raw_data.metabolites})

		# These updates can be dependent on metabolite masses
		self._all_mass.update(self._build_polymerized_subunit_masses(sim_data))
		self._all_mass.update(self._build_rna_masses(raw_data, sim_data))
		self._all_mass.update(self._build_protein_masses(raw_data, sim_data))
		self._all_mass.update(self._build_full_chromosome_mass(raw_data, sim_data))
		self._all_mass.update(
			{x['id']: np.sum(x['mw'])
				for x in itertools.chain(raw_data.proteinComplexes, raw_data.modifiedForms)}
			)

		self._mass_units = units.g / units.mol
		self._location_tag = re.compile(r'\[[a-z]\]')

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
			end_group_mw = self._all_mass[end_group_id[:-3]]
			polymerized_subunit_mws = [
				self._all_mass[met_id[:-3]] - end_group_mw
				for met_id in subunit_ids]

			# Add to dictionary with prefixed ID
			polymerized_subunit_masses.update({
				POLYMERIZED_FRAGMENT_PREFIX + subunit_id[:-3]: mw
					for subunit_id, mw in zip(subunit_ids, polymerized_subunit_mws)
				})

		add_subunit_masses(
			sim_data.moleculeGroups.ntps, sim_data.moleculeIds.ppi)
		add_subunit_masses(
			sim_data.moleculeGroups.dntps, sim_data.moleculeIds.ppi)
		add_subunit_masses(
			sim_data.moleculeGroups.amino_acids, sim_data.moleculeIds.water)

		return polymerized_subunit_masses

	def _build_rna_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights of RNAs keyed with the RNA IDs.
		Molecular weights are calculated from the RNA sequence and the weights
		of polymerized NTPs.
		"""
		# Get RNA nucleotide compositions
		rna_seqs = self.get_sequence([rna['id'] for rna in raw_data.rnas])
		nt_counts = []
		for seq in rna_seqs:
			nt_counts.append(
				[seq.count(letter) for letter in sim_data.ntp_code_to_id_ordered.keys()])
		nt_counts = np.array(nt_counts)

		# Calculate molecular weights
		ppi_mw = self._all_mass[sim_data.moleculeIds.ppi[:-3]]
		polymerized_ntp_mws = np.array([
			self._all_mass[met_id[:-3]] for met_id in sim_data.moleculeGroups.polymerized_ntps
			])

		mws = nt_counts.dot(polymerized_ntp_mws) + ppi_mw  # Add end weight

		return {rna['id']: mw for (rna, mw) in zip(raw_data.rnas, mws)}

	def _build_protein_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights of protein monomers keyed with
		the protein IDs. Molecular weights are calculated from the protein
		sequence and the weights of polymerized amino acids.
		"""
		# Get protein amino acid compositions
		protein_seqs = self.get_sequence(
			[protein['id'] for protein in raw_data.proteins])
		aa_counts = []
		for seq in protein_seqs:
			aa_counts.append(
				[seq.count(letter) for letter in sim_data.amino_acid_code_to_id_ordered.keys()])
		aa_counts = np.array(aa_counts)

		# Calculate molecular weights
		water_mw = self._all_mass[sim_data.moleculeIds.water[:-3]]
		polymerized_aa_mws = np.array(
			[self._all_mass[met_id[:-3]] for met_id in sim_data.moleculeGroups.polymerized_amino_acids]
			)
		mws = aa_counts.dot(polymerized_aa_mws) + water_mw  # Add end weight

		return {protein['id']: mw for (protein, mw) in zip(raw_data.proteins, mws)}

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
			self._all_mass[met_id[:-3]] for met_id in sim_data.moleculeGroups.polymerized_dntps
			])
		mw = np.dot(
			forward_strand_nt_counts + reverse_strand_nt_counts,
			polymerized_dntp_mws)

		return {sim_data.moleculeIds.full_chromosome[:-3]: mw}

	def _build_locations(self, raw_data, sim_data):
		locationDict = {
			item["id"]: list(item["location"])
			for item in itertools.chain(
				raw_data.proteinComplexes,
				raw_data.modifiedForms)}

		# RNAs and full chromosomes only localize to the cytosol
		locationDict.update({
			rna['id']: ['c'] for rna in raw_data.rnas
			})
		locationDict.update({
			sim_data.moleculeIds.full_chromosome[:-3]: ['c']
			})

		# Proteins localize to the single compartment specified in raw data for
		# each protein
		locationDict.update({
			protein['id']: [protein['location']] for protein in raw_data.proteins
			})

		# Metabolites and polymerized subunits can localize to all compartments
		# being modeled
		all_compartments = [comp['abbrev'] for comp in raw_data.compartments]
		locationDict.update({
			met['id']: all_compartments for met in raw_data.metabolites
			})
		locationDict.update({
			subunit_id[:-3]: all_compartments
			for subunit_id in sim_data.moleculeGroups.polymerized_subunits
			})

		# Proteins localize to the single compartment specified in raw data for
		# each protein
		locationDict.update({
			protein['id']: [protein['location']] for protein in raw_data.proteins
			})

		self._locationDict = locationDict
