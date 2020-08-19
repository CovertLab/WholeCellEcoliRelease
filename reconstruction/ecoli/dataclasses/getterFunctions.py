"""
SimulationData getter functions

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import absolute_import, division, print_function

import itertools
import re
from typing import List, Union

import numpy as np

# Unit imports
from wholecell.utils import units


class getterFunctions(object):
	""" getterFunctions """

	def __init__(self, raw_data, sim_data):
		self._build_rna_sequences(raw_data)
		self._build_all_masses(raw_data, sim_data)
		self._build_locations(raw_data)

	def get_rna_sequence(self, ids):
		# type: (Union[List[str], np.ndarray]) -> List[str]
		assert isinstance(ids, (list, np.ndarray))
		return [self._rna_sequences[rna_id] for rna_id in ids]

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
		self._rna_sequences = {}

		for rna in raw_data.rnas:
			gene_index = rna_id_to_gene_index[rna['id']]
			coordinate = coordinate_list[gene_index]

			# Parse genome sequence to get RNA sequence
			if direction_list[gene_index] == '+':
				seq = str(
					genome_sequence[coordinate:coordinate + gene_lengths[gene_index]]
					).replace('T','U')
			else:
				seq = str(
					genome_sequence[coordinate - gene_lengths[gene_index] + 1:coordinate + 1].reverse_complement()
					).replace('T', 'U')

			self._rna_sequences[rna['id']] = seq

	def _build_all_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights keyed with the IDs of molecules.
		Depending on the type of the molecule, weights are either pulled
		directly from raw data, or dynamically calculated from the existing
		data.
		"""
		all_mass = {}

		all_mass.update(self._build_rna_masses(raw_data, sim_data))
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.proteins})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.proteinComplexes})
		all_mass.update({x['id']: np.sum(x['mw7.2']) for x in raw_data.metabolites})
		all_mass.update({x['id']: np.sum(x['mw7.2']) for x in raw_data.modifiedForms})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.polymerized})
		all_mass.update({x['id']: np.sum(x['mw7.2']) for x in raw_data.water})
		all_mass.update({x['id']: np.sum(x['mw']) for x in raw_data.full_chromosome})

		self._all_mass = all_mass
		self._mass_units = units.g / units.mol
		self._location_tag = re.compile(r'\[[a-z]\]')

	def _build_rna_masses(self, raw_data, sim_data):
		"""
		Builds dictionary of molecular weights of RNAs keyed with the RNA IDs.
		"""
		# Get NTP IDs and 1-letter abbreviations
		ntp_ids = sim_data.moleculeGroups.ntpIds
		ntp_abbreviations = [ntp_id[0] for ntp_id in ntp_ids]

		# Get RNA nucleotide compositions
		rna_seqs = self.get_rna_sequence([rna['id'] for rna in raw_data.rnas])
		nt_counts = []
		for seq in rna_seqs:
			nt_counts.append([seq.count(letter) for letter in ntp_abbreviations])
		nt_counts = np.array(nt_counts)

		# Calculate molecular weights
		metabolite_id_to_mw = {
			f'{met["id"]}[c]': met['mw7.2'] for met in raw_data.metabolites
			}
		ppi_mw = metabolite_id_to_mw[sim_data.moleculeIds.ppi]
		polymerized_ntp_mws = np.array(
			[metabolite_id_to_mw[met_id] for met_id in ntp_ids]
			) - ppi_mw

		mws = nt_counts.dot(polymerized_ntp_mws) + ppi_mw  # Add end weight

		return {rna['id']: mw for (rna, mw) in zip(raw_data.rnas, mws)}


	def _build_locations(self, raw_data):
		locationDict = {
			item["id"]: list(item["location"])
			for item in itertools.chain(
				raw_data.proteins,
				raw_data.proteinComplexes,
				raw_data.metabolites,
				raw_data.polymerized,
				raw_data.water,
				raw_data.modifiedForms)}

		# RNAs only localize to the cytosol
		locationDict.update({
			rna['id']: ['c'] for rna in raw_data.rnas
			})

		self._locationDict = locationDict
