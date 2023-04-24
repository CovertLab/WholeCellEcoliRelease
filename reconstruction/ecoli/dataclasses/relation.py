"""
SimulationData relation functions
"""

from __future__ import absolute_import, division, print_function

import warnings
import numpy as np
from Bio.Seq import Seq
from scipy.optimize import minimize
import json
import copy

from wholecell.utils import units
from wholecell.utils.polymerize import polymerize

print_optimization = False

class Relation(object):
	""" Relation """

	def __init__(self, raw_data, sim_data):
		self._build_cistron_to_monomer_mapping(raw_data, sim_data)
		self._build_monomer_to_mRNA_cistron_mapping(raw_data, sim_data)
		self._build_monomer_to_tu_mapping(raw_data, sim_data)
		self._build_RNA_to_tf_mapping(raw_data, sim_data)
		self._build_tf_to_RNA_mapping(raw_data, sim_data)

		# Relate tRNAs, codons, and translation
		self._build_codon_sequences(raw_data, sim_data)
		self._build_codon_based_translation(raw_data, sim_data)
		self._build_codon_dependent_trna_charging(raw_data, sim_data)
		self._build_trna_charging_kinetics(raw_data, sim_data)

	def _build_cistron_to_monomer_mapping(self, raw_data, sim_data):
		"""
		Build a vector that can map vectors that describe a property for RNA
		cistrons into a vector that describes the same property for the
		corresponding monomers if used as an index array. Assumes that each
		monomer maps to a single RNA cistron (A single RNA can map to multiple
		monomers).

		e.g.
		monomer_property = RNA_cistron_property[
			sim_data.relation.cistron_to_monomer_mapping]
		"""
		# Map cistron IDs to indexes given in cistron_data (rnas.tsv)
		cistron_id_to_index = {
			cistron_id: i for i, cistron_id
			in enumerate(sim_data.process.transcription.cistron_data['id'])
			}

		# List the cistron_data indexes of cistron IDs in the order of
		# corresponding cistrons given in monomer_data (proteins.tsv)
		self.cistron_to_monomer_mapping = np.array([
			cistron_id_to_index[cistron_id] for cistron_id
			in sim_data.process.translation.monomer_data['cistron_id']
			])

	def _build_monomer_to_mRNA_cistron_mapping(self, raw_data, sim_data):
		"""
		Builds a sparse matrix that can map vectors that describe a property
		for protein monomers into a vector that describes the same property for
		the corresponding mRNA cistrons if multiplied to the right of the
		original vector. The transformed property must be additive (i.e. if two
		proteins map to the same cistron, the values given for the two proteins
		are added to yield a value for the cistron).

		The full matrix can be returned by calling 
		monomer_to_mRNA_cistron_mapping().
		"""
		# Initialize sparse matrix variables
		self._monomer_to_mRNA_cistron_mapping_i = []
		self._monomer_to_mRNA_cistron_mapping_j = []
		self._monomer_to_mRNA_cistron_mapping_v = []
		self._monomer_to_mRNA_cistron_mapping_shape = (
			len(sim_data.process.translation.monomer_data),
			sim_data.process.transcription.cistron_data['is_mRNA'].sum()
		)

		# Build mapping from mRNA ID to mRNA index
		mRNA_data = sim_data.process.transcription.cistron_data[
			sim_data.process.transcription.cistron_data['is_mRNA']]
		mRNA_id_to_index = {
			mRNA['id']: j for j, mRNA in enumerate(mRNA_data)
		}

		# Build sparse matrix
		for i, monomer in enumerate(sim_data.process.translation.monomer_data):
			self._monomer_to_mRNA_cistron_mapping_i.append(i)
			self._monomer_to_mRNA_cistron_mapping_j.append(
				mRNA_id_to_index[monomer['cistron_id']])
			self._monomer_to_mRNA_cistron_mapping_v.append(1)

	def monomer_to_mRNA_cistron_mapping(self):
		"""
		Returns the full version of the sparse matrix built by
		_build_monomer_to_mRNA_cistron_mapping().

		e.g.
		mRNA_property = sim_data.relation.monomer_to_mRNA_cistron_mapping().T.dot(
			monomer_property)
		"""
		out = np.zeros(self._monomer_to_mRNA_cistron_mapping_shape, dtype=np.float64)
		out[self._monomer_to_mRNA_cistron_mapping_i, self._monomer_to_mRNA_cistron_mapping_j] = self._monomer_to_mRNA_cistron_mapping_v
		return out

	def _build_monomer_to_tu_mapping(self, raw_data, sim_data):
		"""
		Builds a dictionary that maps monomer IDs to a list of all transcription
		unit IDs that the monomer can be translated from.
		"""
		self.monomer_index_to_tu_indexes = {
			i: sim_data.process.transcription.cistron_id_to_rna_indexes(monomer['cistron_id'])
			for i, monomer in enumerate(sim_data.process.translation.monomer_data)
			}

	def _build_RNA_to_tf_mapping(self, raw_data, sim_data):
		"""
		Builds a dictionary that maps RNA IDs to a list of all transcription
		factor IDs that regulate the given RNA. All TFs that target any of the
		constituent cistrons in the RNA are added to each list.
		"""
		cistron_ids = sim_data.process.transcription.cistron_data['id']

		self.rna_id_to_regulating_tfs = {}
		for rna_id in sim_data.process.transcription.rna_data['id']:
			tf_list = []
			for cistron_index in sim_data.process.transcription.rna_id_to_cistron_indexes(rna_id):
				tf_list.extend(
					sim_data.process.transcription_regulation.target_tf.get(
						cistron_ids[cistron_index], []))

			# Remove duplicates and sort
			self.rna_id_to_regulating_tfs[rna_id] = sorted(set(tf_list))

	def _build_tf_to_RNA_mapping(self, raw_data, sim_data):
		"""
		Builds a dictionary that maps transcription factor IDs to a list of all
		RNA IDs that are targeted by the given TF. All RNA transcription units
		that contain any of the cistrons regulated by the TF are added to each
		list.
		"""
		self.tf_id_to_target_RNAs = {}
		for (rna_id, tf_list) in self.rna_id_to_regulating_tfs.items():
			for tf_id in tf_list:
				self.tf_id_to_target_RNAs.setdefault(tf_id, []).append(rna_id)

	def _build_codon_sequences(self, raw_data, sim_data):
		""" Builds the effective codon sequences. """

		def assemble_coding_segments(rna_sequence, coding_segments):
			"""
			Uses coding segments to parse the rna sequence and assemble
			the effective codon sequence.

			rna_sequence (Bio.Seq.Seq): RNA sequence
			coding_segments (List of Lists): describes the start and end
				positions of each coding segment
			"""
			sequence = ''
			for coding_segment in coding_segments:
				start, end = coding_segment
				start -= 1
				end -= 1
				sequence += rna_sequence[start:end+1]
			return sequence

		# Describe coding opal codons (UGA)
		# Note: if more coding UGA codons are discovered in the future,
		# update this dict. Or, optionally retrive this data directly
		# from EcoCyc. At the time of writing this code, it was
		# determined that coding UGA codons are rare and possibly slow
		# to discover instances, so hard-coding was decided.
		coding_opal = {
			'FDNG-MONOMER': 195,
			'FDOG-MONOMER': 195,
			'FORMATEDEHYDROGH-MONOMER': 139,
			}

		# Describe coding segments (describes ribosomal frame shifts and
		# downstream start codons).
		coding_segments = {}
		for rna in raw_data.rnas:

			if len(rna['coding_segments']) == 0:
				continue

			for monomer_id, coding_segment in zip(
					rna['monomer_ids'], rna['coding_segments']):

				if len(coding_segment) == 0:
					continue

				coding_segments[monomer_id] = coding_segment

		# Get sequences
		protein_ids = sim_data.process.translation.monomer_data['id']
		protein_sequences = sim_data.getter.get_sequences(
			[protein_id[:-3] for protein_id in protein_ids])
		rna_sequences = sim_data.getter.get_sequences(
			[rna_id[:-3] for rna_id
			in sim_data.process.transcription.rna_data['id']])

		# Note: Start codons are distinguished from non-starting AUGs.
		self.codons = ([sim_data.molecule_ids.start_codon]
			+ sim_data.molecule_groups.codons)
		codon_to_amino_acid = {codon: [] for codon in self.codons}

		# Build and verify codon sequences
		self._codon_sequences = {}
		for i, protein_id in enumerate(protein_ids):
			protein_id = protein_id[:-3]

			# Get mRNA sequence
			rna_sequence = rna_sequences[self.cistron_to_monomer_mapping[i]]

			# Represent ribosomal frame shifting
			if protein_id in coding_segments:
				rna_sequence = assemble_coding_segments(
					rna_sequence, coding_segments[protein_id])

			# Translate the mRNA sequence
			rna_sequence_translated = rna_sequence.translate()

			# Represent selenocysteine-coding opal codons
			if protein_id in coding_opal:
				opal = coding_opal[protein_id]
				rna_sequence_translated = (rna_sequence_translated[:opal]
					+ 'U' + rna_sequence_translated[opal+1:])

			# Remove stop codons
			stop = rna_sequence_translated.find('*')
			rna_sequence_translated = rna_sequence_translated[:stop]
			rna_sequence = str(rna_sequence[:stop*3])

			# Ensure all sequences begin with Methionine
			# Note: Translation initiation beginning on non-AUG start
			# codons also code for methionine because only initiator
			# tRNAs (charged with formylated methionine) can be
			# recognized by the inititation complex.
			rna_sequence_translated = 'M' + rna_sequence_translated[1:]

			# Get protein sequence
			protein_sequence = protein_sequences[i]

			# Verify codon sequence
			if rna_sequence_translated != protein_sequence:
				warnings.warn('mRNA sequence does not match the protein '
					'sequence for {}'.format(protein_id))
				continue

			# Parse effective mRNA sequence into codons
			codon_sequence = [rna_sequence[j:j+3]
				for j in range(0, len(rna_sequence), 3)]

			# Distinguish start codons
			codon_sequence[0] = sim_data.molecule_ids.start_codon
			self._codon_sequences[protein_id] = codon_sequence

			# Record codon to amino acid interactions
			for codon, amino_acid in zip(codon_sequence, protein_sequence):
				if amino_acid not in codon_to_amino_acid[codon]:
					codon_to_amino_acid[codon].append(amino_acid)

		# Check for overloaded codons
		assert np.all(np.array(
			[len(amino_acids) for amino_acids in codon_to_amino_acid.values()]
			) <= 1)

		# Map from amino acids to codons
		self.amino_acid_to_codons = {amino_acid: []
			for amino_acid in sim_data.molecule_groups.amino_acids}

		for codon, amino_acid in codon_to_amino_acid.items():
			if len(amino_acid) == 0:
				continue

			amino_acid_id = sim_data.amino_acid_code_to_id_ordered[
				amino_acid[0]]

			if codon not in self.amino_acid_to_codons[amino_acid_id]:
				self.amino_acid_to_codons[amino_acid_id].append(codon)

		# Map from codons to amino acids (matrix)
		amino_acids = sim_data.molecule_groups.amino_acids
		self.codons_to_amino_acids = np.zeros(
			(len(amino_acids), len(self.codons)), dtype=np.int8)
		for i, amino_acid in enumerate(amino_acids):
			for codon in self.amino_acid_to_codons[amino_acid]:
				j = self.codons.index(codon)
				self.codons_to_amino_acids[i, j] = 1

	def _build_codon_based_translation(self, raw_data, sim_data):
		"""
		Builds tools used during simulation of polypeptide elongation.

		Notably:
		codon_sequences (analogous to translation_sequences)
		residue_weights_by_codon (analogous to translation_monomer_weights)
		"""

		# Get sequences
		translation = sim_data.process.translation
		protein_ids = translation.monomer_data['id']
		molecule_ids = [protein_id[:-3] for protein_id in protein_ids]
		sequences = [self._codon_sequences[molecule_id]
			for molecule_id in molecule_ids]

		# Calculate the number of columns required to read beyond
		# (towards the C-terminal) the longest mRNA coding region.
		# Note: For polypeptides undergoing N-terminal cleavage of the
		# initial methionine, their sequences (built here) are the
		# pre-cleavage sequence.
		max_len = np.int64(
			# Length of the open reading frame
			+ translation.monomer_data['length'].asNumber().max()

			# Max number of anticipated readings
			+ (translation.max_time_step
				* sim_data.constants.ribosome_elongation_rate_max.asNumber(
					units.aa / units.s)))

		self.codon_sequences = np.full(
			(len(sequences), max_len),
			polymerize.PAD_VALUE, dtype=np.int8)

		# Build mapping
		codon_mapping = {codon:i for i, codon in enumerate(self.codons)}

		# Build translation sequences
		codon_counts = []
		min_length = 1 + max(codon_mapping.values())
		for i, sequence in enumerate(sequences):
			for j, codon in enumerate(sequence):
				self.codon_sequences[i, j] = codon_mapping[codon]

			codon_count = np.bincount(
				self.codon_sequences[i, :j + 1], minlength=min_length)
			codon_counts.append(codon_count)
		self.codon_counts = np.array(codon_counts)

		# Describe residue masses
		residue_weights_by_codon = []
		for col, codon in enumerate(self.codons):
			i = np.where(self.codons_to_amino_acids[:, col])[0][0]
			residue_weights_by_codon.append(translation.translation_monomer_weights[i])
		self.residue_weights_by_codon = np.array(residue_weights_by_codon)

	def _build_codon_dependent_trna_charging(self, raw_data, sim_data):
		"""
		Builds tools used during simulation of tRNA charging.

		Notably:

		"""

		def is_base_pair(codon, anticodon):
			'''
			Returns True if codon and anticodon are base paired, False
			otherwise.

			Note: expects codon and anticodon in 5' to 3' direction.

			Example:
			codon 		5'-ABC-3'
			anticodon 	3'-XYZ-5'

			C and Z are the 'wobble' position
			'''
			for i, (C, A) in enumerate(zip(codon, anticodon[::-1])):

				# Wobble position
				if i == 2:
					return np.any(
						[watson_crick_base_pairing[A] == C,
						wobble_base_pairing[A] == C])

				# Non-wobble positions
				elif watson_crick_base_pairing[A] != C:
					return False

		# Base pairing rules
		watson_crick_base_pairing = {
			'A': 'U',
			'U': 'A',
			'G': 'C',
			'C': 'G'}
		wobble_base_pairing = {
			'A': '',
			'U': 'G',
			'G': 'U',
			'C': ''}

		# Map tRNAs to their anticodons
		rna_data = sim_data.process.transcription.rna_data
		free_trnas = rna_data['id'][rna_data['is_tRNA']]
		anticodons = rna_data['anticodon'][rna_data['is_tRNA']]
		trna_to_anticodon = dict(zip(free_trnas, anticodons))

		# Map free and charged tRNAs
		charged_trnas = sim_data.process.transcription.charged_trna_names
		free_to_charged = dict(zip(free_trnas, charged_trnas))
		self.charged_to_free = dict(zip(charged_trnas, free_trnas))

		# Describe initiator tRNAs
		free_initiators = sim_data.molecule_groups.initiator_trnas
		charged_initiators = [free_to_charged[trna] for trna in free_initiators]

		# Map codons to tRNAs
		transcription = sim_data.process.transcription
		self.codon_to_trnas = {codon: [] for codon in self.codons}
		self.amino_acid_to_trnas = {}

		for i, amino_acid in enumerate(sim_data.molecule_groups.amino_acids):

			# Get trnas
			trnas = free_trnas[np.where(transcription.aa_from_trna[i, :])[0]]
			self.amino_acid_to_trnas[amino_acid] = trnas.tolist()

			# Get codons
			codons = self.amino_acid_to_codons[amino_acid]

			for codon in codons:

				# Start codons are read by initiator trnas
				if codon == sim_data.molecule_ids.start_codon:
					self.codon_to_trnas[codon] = charged_initiators
					continue

				# Elongating AUG codons are read by elongator trnas
				elif codon == 'AUG':
					for trna in trnas:
						if trna not in free_initiators:
							self.codon_to_trnas[codon].append(
								free_to_charged[trna])
					continue

				# ARG trnas argQ, argV, argY, and argZ have inosine (a
				# derivative of adenine) at position 34 of the anticodon
				# stem loop (I34) via post-transcriptional modification,
				# which enables recognition of the cognate codon CGC and
				# wobble codons CGU and CGA.
				elif codon in ['CGA', 'CGC', 'CGU']:
					self.codon_to_trnas[codon] = [
						free_to_charged['argQ-tRNA[c]'],
						free_to_charged['argV-tRNA[c]'],
						free_to_charged['argY-tRNA[c]'],
						free_to_charged['argZ-tRNA[c]'],
						]
					continue

				# ILE trnas ileX and ileY have lysidine (a derivative
				# of cytidine) at position 34 of the anticodon stem loop
				# (L34) via post-transcriptional modification, which
				# enables recognition of the cognate codon AUA.
				elif codon == 'AUA':
					self.codon_to_trnas[codon] = [
						free_to_charged['ileX-tRNA[c]'],
						free_to_charged['RNA0-305[c]'],
						]
					continue

				for trna in trnas:
					if is_base_pair(codon, trna_to_anticodon[trna]):
						self.codon_to_trnas[codon].append(free_to_charged[trna])

				# Assign all tRNAs (that are known to deliver this
				# codon's amino acid) to this codon if no base pairing
				# rule can match this codon to a specific tRNA.
				if len(self.codon_to_trnas[codon]) == 0:
					print(f'Warning: No tRNAs were found to read codon {codon}.'
						f' Assigning all trnas for {amino_acid} to this codon.')
					self.codon_to_trnas[codon] = [free_to_charged[trna]
						for trna in trnas]

		# Map tRNAs to codons
		self.trnas_to_codons = np.zeros(
			(len(self.codons), len(free_trnas)), dtype=np.int8)
		for i, codon in enumerate(self.codons):
			for trna in self.codon_to_trnas[codon]:
				j = np.where(free_trnas == self.charged_to_free[trna])[0][0]
				self.trnas_to_codons[i, j] = 1

		# Describe tRNA-codon pairs
		self.trna_codon_pairs = []
		for col, trna in enumerate(free_trnas):
			for row in np.where(self.trnas_to_codons[:, col])[0]:
				codon = self.codons[row]
				self.trna_codon_pairs.append(f'{trna}_{codon}')

		# Remove lysU; lysU is thought to be used during elevated
		# temperatures, whereas lysS is normally used.
		row = sim_data.molecule_groups.amino_acids.index('LYS[c]')
		col = transcription.synthetase_names.index('LYSU-CPLX[c]')
		self.modified_aa_from_synthetase = np.copy(
			transcription.aa_from_synthetase)
		self.modified_aa_from_synthetase[row, col] = 0

		# Map amino acids to synthetases
		self.amino_acid_to_synthetase = {}
		for i, amino_acid in enumerate(sim_data.molecule_groups.amino_acids):
			synthetase = transcription.synthetase_names[
				np.where(self.modified_aa_from_synthetase[i, :])[0][0]]
			self.amino_acid_to_synthetase[amino_acid] = synthetase

	def _build_trna_charging_kinetics(self, raw_data, sim_data):
		'''
		Parses and stores kinetic parameters recorded in
		trna_charging_kinetics.tsv into the sim_data object.
		'''

		def update(sweep_degree, row):
			synthetase = row['synthetase_id']
			D = self.trna_charging_kinetics[sweep_degree]

			D['synthetase_to_k_cat'].update({synthetase: row['k_cat']})
			D['synthetase_to_K_A'].update({synthetase: row['K_M_amino_acid']})
			D['trna_to_K_T'].update(
				{k: self.conc_unit * v for k,v in row['K_M_trna'].items()})
			D['trna_condition_to_free_fraction'].update(row['f_free'])
			D['synthetase_to_objective'].update({synthetase: row['objective']})

			return

		# Units
		self.conc_unit = units.umol/units.L
		self.rate_unit = self.conc_unit/units.s

		# Store trna charging kinetics
		self.trna_charging_kinetics = {}
		for row in raw_data.optimization.trna_charging_kinetics_solutions:
			synthetase = row['synthetase_id']
			sweep_level = row['sweep_level']

			if sweep_level not in self.trna_charging_kinetics:
				self.trna_charging_kinetics[sweep_level] = {
					'synthetase_to_k_cat': {},
					'synthetase_to_K_A': {},
					'trna_to_K_T': {},
					'trna_condition_to_free_fraction': {},
					'synthetase_to_objective': {},
					}

			if synthetase not in self.trna_charging_kinetics\
				[sweep_level]['synthetase_to_objective']:
				update(sweep_level, row)

			elif row['objective'] < self.trna_charging_kinetics\
					[sweep_level]['synthetase_to_objective'][synthetase]:
				update(sweep_level, row)

		# Store defaults
		default_sweep_level = 4
		self.synthetase_to_k_cat = self.trna_charging_kinetics\
			[default_sweep_level]['synthetase_to_k_cat']
		self.synthetase_to_K_A = self.trna_charging_kinetics\
			[default_sweep_level]['synthetase_to_K_A']
		self.trna_to_K_T = self.trna_charging_kinetics\
			[default_sweep_level]['trna_to_K_T']
		self.trna_condition_to_free_fraction = self.trna_charging_kinetics\
			[default_sweep_level]['trna_condition_to_free_fraction']

		# Store constants
		self.trna_charging_constants = {}
		for row in raw_data.optimization.trna_charging_kinetics_constants:
			key = row['synthetase_id__condition']
			self.trna_charging_constants[key] = {}

			for k, v in row.items():
				if k == key:
					continue
				self.trna_charging_constants[key][k] = v

		# Store tRNA synthetase abundances observed from sample
		# simulations
		self.trna_synthetase_samples = {}
		for row in raw_data.optimization.trna_synthetase_dynamic_range:
			key = row['synthetase_condition']
			self.trna_synthetase_samples[key] = {}
			self.trna_synthetase_samples[key]['mean'] = row['mean']
			self.trna_synthetase_samples[key]['std'] = row['std']
			self.trna_synthetase_samples[key]['min'] = row['min']

		# Store curated kinetics
		self.synthetase_to_max_curated_k_cats = {}
		for row in raw_data.trna_charging_kinetics_curated:
			synthetase = row['Synthetase']
			k_cat_max = row['max k_cat']
			self.synthetase_to_max_curated_k_cats[synthetase] = k_cat_max

		return

	def optimize_trna_charging_kinetics(self, sim_data, cell_specs):
		'''
		Optimizes kinetic parameters describing tRNA charging reactions.
		- Aims to minimize differences between the rates of tRNA
		  aminoacylation and aminoacyl-tRNA utilization (by ribosomes).
		- Holds to the principle that protein content of the cell grows
		  exponentially and doubles at the measured doubling time.
		'''

		def objective(x,
			indexes,
			maps,
			v_codons,
			c_synthetase,
			c_synthetase_min,
			c_amino_acid,
			c_trnas,
			):

			# Pre-conditioning
			x = np.power(10, x)

			# Parse candidate solution
			k_cat = x[indexes['k_cat_index']]
			K_A = x[indexes['K_A_index']]
			K_T = x[indexes['K_T_slice']]
			f = x[indexes['f_slice']]
			min_f = x[indexes['min_f_slice']]

			############################################################

			# Calculate charging rate of each free trna
			saturation_amino_acid = c_amino_acid / (K_A + c_amino_acid)
			relative_trnas = ((maps['f_to_cases'] @ f)
				* c_trnas
				/ (maps['K_T_to_cases'] @ K_T))
			trna_sum = maps['cases_to_trna_sum'] @ relative_trnas
			saturation_trnas = relative_trnas / (1 + trna_sum)
			v_charge = (k_cat
				* c_synthetase
				* saturation_amino_acid
				* saturation_trnas)

			# Calculate distribution of codon reading across trnas
			# Note: columns of codons_to_trnas sum to 1
			c_trnas_charged = (1 - (maps['f_to_cases'] @ f)) * c_trnas
			tile = np.tile(c_trnas_charged, (len(v_codons), 1)).T
			codons_to_trnas = np.where(maps['codon_cases_to_trna_cases'], tile, 0)
			denominator = codons_to_trnas.sum(axis=0)
			denominator[denominator == 0] = 1 # to prevent divide by 0
			codons_to_trnas = np.divide(codons_to_trnas, denominator)

			# Calculate usage rate of each charged trna
			v_usage = codons_to_trnas @ v_codons

			# Steady-state cost
			steady_state_error = np.sum(
				np.square(
					1 - (v_charge / v_usage)
					)
				)

			############################################################

			# Calculate charging rate of each free trna
			relative_trnas = ((maps['f_to_cases'] @ min_f)
				* c_trnas
				/ (maps['K_T_to_cases'] @ K_T))
			trna_sum = maps['cases_to_trna_sum'] @ relative_trnas
			saturation_trnas = relative_trnas / (1 + trna_sum)
			v_charge_min = (k_cat
				* c_synthetase_min
				* saturation_amino_acid
				* saturation_trnas)

			# Calculate distribution of codon reading across trnas
			# Note: columns of codons_to_trnas sum to 1
			c_trnas_charged = (1 - (maps['f_to_cases'] @ min_f)) * c_trnas
			tile = np.tile(c_trnas_charged, (len(v_codons), 1)).T
			codons_to_trnas = np.where(maps['codon_cases_to_trna_cases'], tile, 0)
			denominator = codons_to_trnas.sum(axis=0)
			denominator[denominator == 0] = 1 # to prevent divide by 0
			codons_to_trnas = np.divide(codons_to_trnas, denominator)

			# Calculate usage rate of each charged trna
			v_usage_min = codons_to_trnas @ v_codons

			# Steady-state cost
			min_steady_state_error = np.sum(
				np.square(
					1 - (v_charge_min / v_usage_min)
					)
				)

			############################################################

			# Bounds penalty
			bounds_penalty = sum((1 / (f - 0.05)) + (1 / (0.95 - f)))
			bounds_penalty += sum((1 / (min_f - 0.05)) + (1 / (0.95 - min_f)))

			############################################################

			errors = [

				# Steady state errors
				steady_state_error,
				min_steady_state_error,

				# Penalties
				(w_r * sum(K_T)),
				(w_b * bounds_penalty),

				]

			error = sum(errors)

			return error

		def full_saturation_amino_acid(c_amino_acid, c_synthetase,
				v_codons_total):
			'''
			V_rich / V_poor = ((E_rich / E_poor)
				* (sat_A_rich / sat_A_poor)
				* (sat_T_rich / sat_T_poor))

			Assuming:
			1) sat_A_rich = 1
			2) sat_T_poor = sat_T_rich

			V_rich / V_poor = (E_rich / E_poor) * (1 / sat_A_poor)
			sat_A_poor = (E_rich / E_poor) / (V_rich / V_poor)
			           = A_poor / (K_A + A_poor)
			(K_A + A_poor) = A_poor / (E_rich / E_poor) * (V_rich / V_poor)
			K_A = A_poor * ((V_rich / V_poor) / (E_rich / E_poor)  - 1)
			'''
			A_poor, A_rich = c_amino_acid[condition_start_indexes]
			E_poor, E_rich = c_synthetase[condition_start_indexes]
			V_poor, V_rich = v_codons_total
			K_A = A_poor * ((V_rich / V_poor) / (E_rich / E_poor)  - 1)
			return np.log10(K_A)

		def get_random_initial_solution(c_synthetase, c_amino_acid,
				c_trnas, v_codons_total, maps, n_K_T, n_f, indexes,
				condition_start_indexes, K_T_range,
				assume_full_saturation):
			'''
			Random initialization
			'''

			initial_solution = np.zeros(indexes['n_parameters'], np.float64)

			# K_A
			if assume_full_saturation:
				K_A = full_saturation_amino_acid(
					c_amino_acid, c_synthetase, v_codons_total)

			else:
				K_A = np.random.uniform(
					np.min(np.log10(c_amino_acid)),
					np.max(np.log10(c_amino_acid)),
					1)[0]
			initial_solution[indexes['K_A_index']] = K_A

			# K_T
			K_T = np.random.uniform(K_T_range[0], K_T_range[1], n_K_T)
			initial_solution[indexes['K_T_slice']] = K_T

			# f
			n_f = n_conditions * n_K_T
			f = np.log10(np.random.uniform(0.1, 0.9, n_f))
			# f = np.log10(np.random.uniform(0.1, 0.15, n_f))
			initial_solution[indexes['f_slice']] = f

			# Solve for k_cat
			K_A = np.power(10, K_A)
			K_T = np.power(10, K_T)
			f = np.power(10, f)

			saturation_amino_acid = c_amino_acid / (K_A + c_amino_acid)
			relative_trnas = ((maps['f_to_cases'] @ f)
				* c_trnas
				/ (maps['K_T_to_cases'] @ K_T))
			trna_sum = maps['cases_to_trna_sum'] @ relative_trnas
			saturation_trnas = (trna_sum[condition_start_indexes]
				/ (1 + trna_sum[condition_start_indexes]))
			k_cat = np.max(v_codons_total
				/ c_synthetase[condition_start_indexes]
				/ saturation_amino_acid[condition_start_indexes]
				/ saturation_trnas)
			initial_solution[indexes['k_cat_index']] = np.log10(k_cat)

			# f (min)
			min_f = np.log10(np.random.uniform(0.1, 0.9, n_f))
			# min_f = np.log10(np.random.uniform(0.85, 0.9, n_f))
			initial_solution[indexes['min_f_slice']] = min_f

			return initial_solution

		################################################################
		# Optimize tRNA Charging Kinetics and save results

		trna_charging_kinetics_solutions = [[
			'"synthetase_id"',
			'"sweep_level"',
			'"iteration"',
			'"objective"',
			'"k_cat (1/units.s)"',
			'"K_M_amino_acid (units.umol/units.L)"',
			'"K_M_trna"',
			'"f_free"',
			'"f_free_at_min"',
			]]

		trna_charging_kinetics_constants = [[
			'"synthetase_id__condition"',
			'"synthetase (units.umol/units.L)"',
			'"amino_acid (units.umol/units.L)"',
			'"trnas"',
			'"codons"',
			]]

		# Describe weights
		w_b = 1e-9 # bounds
		w_r = 1e-9 # regularization

		# Bounds
		f_lower_bound = 0.051
		f_upper_bound = 0.949

		# Iterations
		iterations = 100
		viable_solutions = 10

		# Other parameters used during optimization
		K_T_range = [1, 3]
		n_conditions = len(self.conditions)
		n_increments = 5
		sweep_levels = np.arange(2, n_increments + 1)
		upper_bound = 10 # log

		# Optimize each amino acid system individually
		for amino_acid in sim_data.molecule_groups.amino_acids:
			if amino_acid == 'L-SELENOCYSTEINE[c]':
				continue

			# Get molecules
			synthetase = self.amino_acid_to_synthetase[amino_acid]
			trnas = self.amino_acid_to_trnas[amino_acid]
			codons = self.amino_acid_to_codons[amino_acid]

			# Get constants
			(c_synthetase, c_amino_acid, c_trnas, v_codons, v_codons_total,
				trna_charging_kinetics_constants) = self.get_constants(
					synthetase,
					amino_acid,
					trnas,
					codons,
					sim_data,
					cell_specs,
					trna_charging_kinetics_constants,
					)

			# Build maps
			n_K_T, K_T_indexes, codons_to_trnas = self.assign_K_T(trnas, codons)

			cases = []
			for condition in self.conditions:
				for trna in trnas:
					cases.append(f'{trna}__{condition}')
			n_cases = len(cases)

			f_to_cases = np.zeros(
				(n_cases, n_K_T * n_conditions), dtype=np.int64)
			for row, case in enumerate(cases):
				trna, condition = case.split('__')
				K_T_index = K_T_indexes[trnas.index(trna)]
				col = (self.conditions.index(condition) * n_K_T) + K_T_index
				f_to_cases[row, col] = 1

			K_T_to_cases = np.zeros((n_cases, n_K_T), dtype=np.int64)
			for row, case in enumerate(cases):
				trna, condition = case.split('__')
				col = K_T_indexes[trnas.index(trna)]
				K_T_to_cases[row, col] = 1

			cases_to_trna_sum = np.zeros((n_cases, n_cases), dtype=np.int64)
			for row, case in enumerate(cases):
				trna, condition = case.split('__')
				cols = [condition in case for case in cases]
				cases_to_trna_sum[row, cols] = 1

			codon_cases_to_trna_cases = np.zeros(
				(n_cases, len(codons) * n_conditions), dtype=np.bool)
			for i in range(n_conditions):
				rows = slice(i * len(trnas), (i+1) * len(trnas))
				cols = slice(i * len(codons), (i+1) * len(codons))
				codon_cases_to_trna_cases[rows, cols] = codons_to_trnas

			maps = {
				'f_to_cases': f_to_cases,
				'K_T_to_cases': K_T_to_cases,
				'cases_to_trna_sum': cases_to_trna_sum,
				'codon_cases_to_trna_cases': codon_cases_to_trna_cases,
				}

			condition_start_indexes = np.arange(0, n_cases, len(trnas))
			condition_indexes = {k: np.arange(0, n_K_T) + (k * n_K_T)
				for k in range(n_conditions)}

			# Set indexes
			n_f = n_conditions * n_K_T
			k_cat_index = 0
			K_A_index = 1
			K_T_slice = slice(2, 2 + n_K_T)
			f_slice = slice(K_T_slice.stop, K_T_slice.stop + n_f)
			min_f_slice = slice(f_slice.stop, f_slice.stop + n_f)
			n_parameters = min_f_slice.stop

			indexes = {
				'k_cat_index': k_cat_index,
				'K_A_index': K_A_index,
				'K_T_slice': K_T_slice,
				'f_slice': f_slice,
				'min_f_slice': min_f_slice,
				'condition_indexes': condition_indexes,
				'n_parameters': n_parameters,
				}

			# Bounds
			K_A = full_saturation_amino_acid(
				c_amino_acid, c_synthetase, v_codons_total)

			bounds_1 = [() for n in range(n_parameters)]
			bounds_1[k_cat_index] = (-2, upper_bound)
			bounds_1[K_A_index] = (K_A, K_A)
			bounds_1[K_T_slice] = [(-2, upper_bound)
				for n in range(n_K_T)]

			bounds_1[f_slice] = [(
				np.log10(f_lower_bound),
				np.log10(f_upper_bound),
				) for n in range(n_f)]
			bounds_1[min_f_slice] = [(
				np.log10(f_lower_bound),
				np.log10(f_upper_bound),
				) for n in range(n_f)]

			bounds_2 = bounds_1[:]
			bounds_2[K_A_index] = (-2, upper_bound)

			# Get the minimum trna synthetase concentration observed
			# from sample simulations
			synthetase_mins = []
			for condition in self.conditions:
				key = f'{synthetase}__{condition}'
				synthetase_mins.append(
					self.trna_synthetase_samples[key]['min']
					.asNumber(self.conc_unit)
					)
			synthetase_mins = np.array(synthetase_mins)

			# Sweep over the minimum tRNA synthetase value
			for sweep_level in sweep_levels:
				if print_optimization:
					print('='*20)
					print(amino_acid, sweep_level)
					print('='*20)

				# Calculate the minimum trna synthetase concentration
				# used in this sweep level
				c_synthetase_min = []
				for estimated_minimum in (sweep_level
						/ n_increments
						* synthetase_mins):
					for i in range(len(trnas)):
						c_synthetase_min.append(estimated_minimum)
				c_synthetase_min = np.array(c_synthetase_min)

				# Determine whether full amino acid saturation holds
				initial_solution = get_random_initial_solution(
					c_synthetase, c_amino_acid, c_trnas, v_codons_total,
					maps, n_K_T, n_f, indexes, condition_start_indexes,
					K_T_range, assume_full_saturation=True)

				result = minimize(
					objective,
					initial_solution,
					method='Powell',
					bounds=bounds_1,
					args=(
						indexes,
						maps,
						v_codons,
						c_synthetase,
						c_synthetase_min,
						c_amino_acid,
						c_trnas,
						),
					options={
						'maxiter': 1000,
						'ftol': 1e-5,
						},
					)

				x = np.power(10, result.x)
				K_A = x[K_A_index]
				saturation_amino_acid = c_amino_acid / (K_A + c_amino_acid)
				assume_full_saturation = np.isclose(
					saturation_amino_acid[condition_start_indexes][1],
					0.99, atol=1e-2)

				if print_optimization:
					print(f'Assuming full AA saturation: {assume_full_saturation}')

				# Solve
				min_objective = np.inf
				viable_solution = 0
				iteration = 0
				while (iteration < iterations) or (viable_solution < viable_solutions):

					# Option 1: Assume full saturation of amino acids
					if assume_full_saturation:

						# Initialize random initial solution
						initial_solution = get_random_initial_solution(
							c_synthetase, c_amino_acid, c_trnas, v_codons_total,
							maps, n_K_T, n_f, indexes, condition_start_indexes,
							K_T_range, assume_full_saturation=True)

						# Optimize
						result = minimize(
							objective,
							initial_solution,
							method='Powell',
							bounds=bounds_1,
							args=(
								indexes,
								maps,
								v_codons,
								c_synthetase,
								c_synthetase_min,
								c_amino_acid,
								c_trnas,
								),
							options={
								'maxiter': 1000,
								'ftol': 1e-5,
								},
							)


					# Option 2: Don't assume full saturation of amino acids
					else:

						# Initialize random initial solution
						initial_solution = get_random_initial_solution(
							c_synthetase, c_amino_acid, c_trnas, v_codons_total,
							maps, n_K_T, n_f, indexes, condition_start_indexes,
							K_T_range, assume_full_saturation=False)

						# Optimize
						result = minimize(
							objective,
							initial_solution,
							method='Powell',
							bounds=bounds_2,
							args=(
								indexes,
								maps,
								v_codons,
								c_synthetase,
								c_synthetase_min,
								c_amino_acid,
								c_trnas,
								),
							options={
								'maxiter': 1000,
								'ftol': 1e-5,
								},
							)

					# Retrieve solution
					x = np.power(10, result.x)

					k_cat = x[k_cat_index]
					K_A = x[K_A_index]
					K_T = x[K_T_slice]
					f = x[f_slice]
					min_f = x[min_f_slice]

					# Calculate charging rate of each free trna
					saturation_amino_acid = c_amino_acid / (K_A + c_amino_acid)
					relative_trnas = ((maps['f_to_cases'] @ f)
						* c_trnas
						/ (maps['K_T_to_cases'] @ K_T))
					trna_sum = maps['cases_to_trna_sum'] @ relative_trnas
					saturation_trnas = relative_trnas / (1 + trna_sum)
					v_charge = (k_cat
						* c_synthetase
						* saturation_amino_acid
						* saturation_trnas)
					sat_T = trna_sum / (1 + trna_sum)

					# Calculate distribution of codon reading across trnas
					# Note: columns of codons_to_trnas sum to 1
					c_trnas_charged = (1 - (maps['f_to_cases'] @ f)) * c_trnas
					tile = np.tile(c_trnas_charged, (len(v_codons), 1)).T
					codons_to_trnas = np.where(maps['codon_cases_to_trna_cases'], tile, 0)
					denominator = codons_to_trnas.sum(axis=0)
					denominator[denominator == 0] = 1 # to prevent divide by 0
					codons_to_trnas = np.divide(codons_to_trnas, denominator)

					# Calculate usage rate of each charged trna
					v_usage = codons_to_trnas @ v_codons

					# Calculate charging rate of each free trna
					relative_trnas = ((maps['f_to_cases'] @ min_f)
						* c_trnas
						/ (maps['K_T_to_cases'] @ K_T))
					trna_sum = maps['cases_to_trna_sum'] @ relative_trnas
					saturation_trnas = relative_trnas / (1 + trna_sum)
					v_charge_min = (k_cat
						* c_synthetase_min
						* saturation_amino_acid
						* saturation_trnas)
					sat_T_min = trna_sum / (1 + trna_sum)

					# Calculate distribution of codon reading across trnas
					# Note: columns of codons_to_trnas sum to 1
					c_trnas_charged = (1 - (maps['f_to_cases'] @ min_f)) * c_trnas
					tile = np.tile(c_trnas_charged, (len(v_codons), 1)).T
					codons_to_trnas = np.where(maps['codon_cases_to_trna_cases'], tile, 0)
					denominator = codons_to_trnas.sum(axis=0)
					denominator[denominator == 0] = 1 # to prevent divide by 0
					codons_to_trnas = np.divide(codons_to_trnas, denominator)

					# Calculate usage rate of each charged trna
					v_usage_min = codons_to_trnas @ v_codons

					iteration += 1

					# Only save solutions that can maintain steady state
					if (np.max(np.abs(v_charge - v_usage)) > 1e-3
						or np.isnan(np.max(np.abs(v_charge - v_usage)))
						or np.max(np.abs(v_charge_min - v_usage_min)) > 1e-3
						or np.isnan(np.max(np.abs(v_charge_min - v_usage_min)))
						):
						continue

					if result.fun < min_objective:
						min_objective = result.fun

						if print_optimization:
							print('='*26)
							print('\nobj\t', f'{result.fun:.1e}')
							print('-'*26)
							print('k_cat\t', f'{k_cat:.2f}')
							print('K_A\t', f'{K_A:.2f}')
							print('K_T\t', np.round(K_T, 2))
							print('sat A\t',
								f'{saturation_amino_acid[0]:.2f}',
								f'\t{saturation_amino_acid[-1]:.2f}')
							print('sat T\t',
								f'{sat_T[0]:.2f}',
								f'\t{sat_T[-1]:.2f}')

							print('-'*8, 'at exp E', '-'*8)
							print('fc\t', np.round(1-f, 2))
							print('v\t',
								f'{v_charge[:len(trnas)].sum():.2f}',
								f'{v_codons[:len(codons)].sum():.2f}')

							print('-'*8, 'at min E', '-'*8)
							print('fc\t', np.round(1-min_f, 2))
							print('v\t',
								f'{v_charge_min[:len(trnas)].sum():.2f}',
								f'{v_codons[:len(codons)].sum():.2f}')

							print('v diff\t\t', f'{np.max(np.abs(v_charge - v_usage)):.2e}')
							print('v diff min\t', f'{np.max(np.abs(v_charge_min - v_usage_min)):.2e}')

					# Only save solutions that do not have nans
					if (np.isnan(result.fun)
							or np.isinf(result.fun)
							or np.any(np.isnan(np.power(10, result.x)))
							or np.any(np.isinf(np.power(10, result.x)))):
						continue

					trna_charging_kinetics_solutions.append([
						f'"{synthetase}"',
						sweep_level,
						iteration,
						result.fun,
						k_cat,
						K_A,
						json.dumps(dict(zip(trnas, K_T[K_T_indexes]))),
						json.dumps(dict(zip(cases, maps['f_to_cases'] @ f))),
						json.dumps(dict(zip(cases, maps['f_to_cases'] @ min_f))),
						])
					viable_solution += 1

		return trna_charging_kinetics_solutions, trna_charging_kinetics_constants

	def get_constants(self, synthetase, amino_acid, trnas, codons,
			sim_data, cell_specs, trna_charging_kinetics_constants,
			select_condition=None):
		'''
		Retrives constants values of the objective function, which are:
		- concentration of trna synthetase enzymes (uM)
		- concentration of amino acids (uM)
		- concentration of total trnas (uM)
		- rate of codon reading (uM/s)
		'''
		c_synthetase_array = []
		c_amino_acid_array = []
		c_trnas_array = []
		v_codons_array = []
		v_codons_total = []

		i_codons = [self.codons.index(codon) for codon in codons]

		for condition in self.conditions:

			# Get amino acid targets used in metabolism
			media_id = sim_data.conditions[condition]['nutrients']
			c_amino_acid = (sim_data
				.process
				.metabolism
				.concentration_updates
				.concentrations_based_on_nutrients(media_id=media_id)
				[amino_acid]).asNumber(self.conc_unit)

			# Get codon reading rate
			v_codons = (sim_data.codon_read_rate[media_id][i_codons]
				).asNumber(self.rate_unit)
			v_codons_array += v_codons.tolist()
			v_codons_total.append(v_codons.sum())

			# Get tRNA synthetase concentration observed from sample
			# simulations
			c_synthetase = (
				self.trna_synthetase_samples[f'{synthetase}__{condition}']['mean']
				).asNumber(self.conc_unit)

			# Get tRNA abundances from average container
			volume = (cell_specs[condition]['avgCellDryMassInit']
				/ sim_data.mass.cell_dry_mass_fraction
				/ sim_data.constants.cell_density)
			to_conc = 1 / sim_data.constants.n_avogadro / volume
			container = cell_specs[condition]['bulkAverageContainer']

			trna_to_conc = {}
			for trna in trnas:
				c_trna = (to_conc
					* container.count(trna)
					).asNumber(self.conc_unit)

				trna_to_conc[trna] = c_trna
				c_trnas_array.append(c_trna)
				c_synthetase_array.append(c_synthetase)
				c_amino_acid_array.append(c_amino_acid)

			# Record
			trna_charging_kinetics_constants.append([
				f'"{synthetase}__{condition}"',
				c_synthetase,
				c_amino_acid,
				json.dumps(trna_to_conc),
				json.dumps(dict(zip(codons, v_codons))),
				])

		c_synthetase_array = np.array(c_synthetase_array)
		c_amino_acid_array = np.array(c_amino_acid_array)
		c_trnas_array = np.array(c_trnas_array)
		v_codons_array = np.array(v_codons_array)
		v_codons_total = np.array(v_codons_total)

		return (c_synthetase_array,
			c_amino_acid_array,
			c_trnas_array,
			v_codons_array,
			v_codons_total,
			trna_charging_kinetics_constants)

	def assign_K_T(self, trnas, codons):
		'''
		Assign K_Ts (the Michaelis-Menten constant describe the affinity
		between tRNA synthetase enzymes and their tRNA substrates) tRNAs

		tRNAs that read the same set of codons are assigned the same K_T
		'''
		# Group tRNAs by the codon(s) they share
		trna_groups = list(set([tuple(self.codon_to_trnas[codon])
			for codon in codons]))

		# Sort by group size
		order = np.argsort([len(group) for group in trna_groups])
		trna_groups  = np.array(trna_groups, dtype=object)[order]

		# Assign unique K_T's according to primary group affiliation
		K_T_indexes = -1 * np.ones(len(trnas), dtype=np.int8)
		for i, group in enumerate(trna_groups):
			for trna in group:
				trna = self.charged_to_free[trna]
				if K_T_indexes[trnas.index(trna)] == -1:
					K_T_indexes[trnas.index(trna)] = i
		n_K_T = max(K_T_indexes) + 1

		# Map codons to tRNAs
		codons_to_trnas = np.zeros((len(trnas), len(codons)), dtype=np.int64)
		for col, codon in enumerate(codons):
			for trna in self.codon_to_trnas[codon]:
				row = trnas.index(self.charged_to_free[trna])
				codons_to_trnas[row, col] = 1

		return n_K_T, K_T_indexes, codons_to_trnas
