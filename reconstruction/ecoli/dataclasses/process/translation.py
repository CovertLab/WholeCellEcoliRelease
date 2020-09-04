"""
SimulationData for translation process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/09/2015
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import six

from wholecell.sim.simulation import MAX_TIME_STEP
from wholecell.utils import data, units
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import make_elongation_rates


PROCESS_MAX_TIME_STEP = 2.


class Translation(object):
	""" Translation """

	def __init__(self, raw_data, sim_data):
		self.max_time_step = min(MAX_TIME_STEP, PROCESS_MAX_TIME_STEP)

		self._buildMonomerData(raw_data, sim_data)
		self._buildTranslation(raw_data, sim_data)
		self._buildTranslationEfficiency(raw_data, sim_data)
		self._build_elongation_rates(raw_data, sim_data)

	def __getstate__(self):
		"""Return the state to pickle with translationSequences removed and
		only storing data from translationSequences with pad values stripped.
		"""

		state = data.dissoc_strict(self.__dict__, ('translationSequences',))
		state['sequences'] = np.array([
			seq[seq != polymerize.PAD_VALUE]
			for seq in self.translationSequences], dtype=object)
		state['sequence_shape'] = self.translationSequences.shape
		return state

	def __setstate__(self, state):
		"""Restore translationSequences and remove processed versions of the data."""
		sequences = state.pop('sequences')
		sequence_shape = state.pop('sequence_shape')
		self.__dict__.update(state)

		self.translationSequences = np.full(sequence_shape, polymerize.PAD_VALUE, dtype=np.int8)
		for i, seq in enumerate(sequences):
			self.translationSequences[i, :len(seq)] = seq

	def _buildMonomerData(self, raw_data, sim_data):
		# Get protein IDs with compartments
		protein_ids = [protein['id'] for protein in raw_data.proteins]
		protein_compartments = sim_data.getter.getLocation(protein_ids)
		assert all([len(loc) == 1 for loc in protein_compartments])
		protein_ids_with_compartments = [
			f'{protein_id}[{loc[0]}]' for (protein_id, loc)
			in zip(protein_ids, protein_compartments)
			]
		n_proteins = len(protein_ids)

		# Get mappings from monomer IDs to RNA IDs
		monomer_id_to_rna_id = {
			rna['monomerId']: rna['id'] for rna in raw_data.rnas}

		# Get RNA IDs with compartments
		rna_ids = [
			monomer_id_to_rna_id[protein['id']]
			for protein in raw_data.proteins]
		rna_compartments = sim_data.getter.getLocation(rna_ids)
		assert all([len(loc) == 1 for loc in rna_compartments])
		rna_ids_with_compartments = [
			f'{rna_id}[{loc[0]}]'
			for (rna_id, loc) in zip(rna_ids, rna_compartments)]

		# Get lengths and amino acids counts of each protein
		protein_seqs = sim_data.getter.get_sequence(protein_ids)
		lengths = [len(seq) for seq in protein_seqs]
		aa_counts = [
			[seq.count(aa) for aa in sim_data.amino_acid_code_to_id_ordered.keys()]
			for seq in protein_seqs]
		n_amino_acids = len(sim_data.amino_acid_code_to_id_ordered)

		# Get molecular weights
		mws = sim_data.getter.getMass(protein_ids).asNumber(units.g/units.mol)

		# Calculate degradation rates based on N-rule
		deg_rate_units = 1 / units.s
		n_end_rule_deg_rates = {
			row['aa_code']: (np.log(2)/(row['half life'])).asNumber(deg_rate_units)
			for row in raw_data.protein_half_lives_n_end_rule}
		slow_deg_rate = min(n_end_rule_deg_rates.values())

		# Build list of ribosomal proteins
		# Give all ribosomal proteins the slowAA rule
		ribosomalProteins = []
		ribosomalProteins.extend([x[:-3] for x in sim_data.moleculeGroups.s30_proteins])
		ribosomalProteins.extend([x[:-3] for x in sim_data.moleculeGroups.s50_proteins])

		# Get degradation rates from measured protein half lives
		measured_deg_rates = {
			p['id']: (np.log(2) / p['half life']).asNumber(deg_rate_units)
			for p in raw_data.protein_half_lives_measured
			}

		degRate = np.zeros(len(raw_data.proteins))
		for i, protein in enumerate(raw_data.proteins):
			if protein['id'] in measured_deg_rates:
				degRate[i] = measured_deg_rates[protein['id']]
			elif protein['id'] not in ribosomalProteins:
				degRate[i] = n_end_rule_deg_rates[protein['seq'][0]]
			else:
				degRate[i] = slow_deg_rate

		max_protein_id_length = max(
			len(protein_id) for protein_id in protein_ids_with_compartments)
		max_rna_id_length = max(
			len(rna_id) for rna_id in rna_ids_with_compartments)
		monomerData = np.zeros(
			n_proteins,
			dtype = [
				('id', 'U{}'.format(max_protein_id_length)),
				('rnaId', 'U{}'.format(max_rna_id_length)),
				('degRate', 'f8'),
				('length', 'i8'),
				('aaCounts', '{}i8'.format(n_amino_acids)),
				('mw', 'f8'),
				]
			)

		monomerData['id'] = protein_ids_with_compartments
		monomerData['rnaId'] = rna_ids_with_compartments
		monomerData['degRate'] = degRate
		monomerData['length'] = lengths
		monomerData['aaCounts'] = aa_counts
		monomerData['mw'] = mws

		field_units = {
			'id'		:	None,
			'rnaId'		:	None,
			'degRate'	:	deg_rate_units,
			'length'	:	units.aa,
			'aaCounts'	:	units.aa,
			'mw'		:	units.g / units.mol,
			}

		self.monomerData = UnitStructArray(monomerData, field_units)
		self.n_monomers = len(self.monomerData)

	def _buildTranslation(self, raw_data, sim_data):
		sequences = sim_data.getter.get_sequence(
			[protein['id'] for protein in raw_data.proteins])

		max_len = np.int64(
			self.monomerData["length"].asNumber().max()
			+ self.max_time_step * sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
			)

		self.translationSequences = np.full((len(sequences), max_len), polymerize.PAD_VALUE, dtype=np.int8)
		aa_ids_single_letter = six.viewkeys(sim_data.amino_acid_code_to_id_ordered)
		aaMapping = {aa: i for i, aa in enumerate(aa_ids_single_letter)}
		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.translationSequences[i, j] = aaMapping[letter]

		aaIDs = list(sim_data.amino_acid_code_to_id_ordered.values())

		self.translationMonomerWeights = (
			(
				sim_data.getter.getMass(aaIDs)
				- sim_data.getter.getMass([sim_data.moleculeIds.water])
				)
			/ sim_data.constants.nAvogadro
			).asNumber(units.fg)
		self.translationEndWeight = (sim_data.getter.getMass([sim_data.moleculeIds.water]) / sim_data.constants.nAvogadro).asNumber(units.fg)

	def _buildTranslationEfficiency(self, raw_data, sim_data):
		monomer_ids = [protein["id"] for protein in raw_data.proteins]

		# Get mappings from monomer IDs to gene IDs
		monomer_id_to_rna_id = {
			rna['monomerId']: rna['id'] for rna in raw_data.rnas}
		rna_id_to_gene_id = {
			gene['rnaId']: gene['id'] for gene in raw_data.genes}
		monomer_id_to_gene_id = {
			monomer_id: rna_id_to_gene_id[monomer_id_to_rna_id[monomer_id]]
			for monomer_id in monomer_ids}

		# Get mappings from gene IDs to translation efficiencies
		gene_id_to_trl_eff = {
			x["geneId"]: x["translationEfficiency"]
			for x in raw_data.translationEfficiency
			if type(x["translationEfficiency"]) == float}

		trl_effs = []
		for monomer_id in monomer_ids:
			gene_id = monomer_id_to_gene_id[monomer_id]

			if gene_id in gene_id_to_trl_eff:
				trl_effs.append(gene_id_to_trl_eff[gene_id])
			else:
				trl_effs.append(np.nan)

		# If efficiency is unavailable, the average of existing effciencies
		# is used
		self.translationEfficienciesByMonomer = np.array(trl_effs)
		self.translationEfficienciesByMonomer[
			np.isnan(self.translationEfficienciesByMonomer)
			] = np.nanmean(self.translationEfficienciesByMonomer)

	def _build_elongation_rates(self, raw_data, sim_data):
		protein_ids = self.monomerData['id']
		ribosomal_protein_ids = sim_data.moleculeGroups.rProteins

		protein_indexes = {
			protein: index
			for index, protein in enumerate(protein_ids)}

		ribosomal_proteins = {
			rprotein: protein_indexes.get(rprotein, -1)
			for rprotein in ribosomal_protein_ids}

		self.rprotein_indexes = np.array([
			index
			for index in sorted(ribosomal_proteins.values())
			if index >= 0], dtype=np.int64)

		self.basal_elongation_rate = sim_data.constants.ribosomeElongationRateBasal.asNumber(units.aa / units.s)
		self.max_elongation_rate = sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
		self.elongation_rates = np.full(
			self.n_monomers,
			self.basal_elongation_rate,
			dtype=np.int64)

		self.elongation_rates[self.rprotein_indexes] = self.max_elongation_rate

	def make_elongation_rates(
			self,
			random,
			base,
			time_step,
			variable_elongation=False):

		return make_elongation_rates(
			random,
			self.n_monomers,
			base,
			self.rprotein_indexes,
			self.max_elongation_rate,
			time_step,
			variable_elongation)
