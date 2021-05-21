"""
SimulationData for transcription process

TODO: add mapping of tRNA to charged tRNA if allowing more than one modified form of tRNA and separate mappings for tRNA and charged tRNA to AA
TODO: handle ppGpp and DksA-ppGpp regulation separately
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import scipy
from scipy import interpolate
import sympy as sp
from typing import cast

from wholecell.sim.simulation import MAX_TIME_STEP
from wholecell.utils import data, units
from wholecell.utils.fitting import normalize
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import make_elongation_rates
from six.moves import zip


PROCESS_MAX_TIME_STEP = 2.
RNA_SEQ_ANALYSIS = "rsem_tpm"
KCAT_ENDO_RNASE = 0.001
ESTIMATE_ENDO_RNASES = 5000
PPGPP_CONC_UNITS = units.umol / units.L
PRINT_VALUES = False  # print values for supplemental table if True


class Transcription(object):
	"""
	SimulationData for the transcription process
	"""

	def __init__(self, raw_data, sim_data):
		self.max_time_step = min(MAX_TIME_STEP, PROCESS_MAX_TIME_STEP)

		self._build_ppgpp_regulation(raw_data, sim_data)
		self._build_rna_data(raw_data, sim_data)
		self._build_transcription(raw_data, sim_data)
		self._build_charged_trna(raw_data, sim_data)
		self._build_attenuation(raw_data, sim_data)
		self._build_elongation_rates(raw_data, sim_data)

	def __getstate__(self):
		"""Return the state to pickle with transcriptionSequences removed and
		only storing data from transcriptionSequences with pad values stripped.
		"""

		state = data.dissoc_strict(self.__dict__, ('transcription_sequences',))
		state['sequences'] = np.array([
			seq[seq != polymerize.PAD_VALUE]
			for seq in self.transcription_sequences], dtype=object)
		state['sequence_shape'] = self.transcription_sequences.shape
		return state

	def __setstate__(self, state):
		"""Restore transcriptionSequences and remove processed versions of the data."""
		sequences = state.pop('sequences')
		sequence_shape = state.pop('sequence_shape')
		self.__dict__.update(state)

		self.transcription_sequences = np.full(sequence_shape, polymerize.PAD_VALUE, dtype=np.int8)
		for i, seq in enumerate(sequences):
			self.transcription_sequences[i, :len(seq)] = seq

	def _build_ppgpp_regulation(self, raw_data, sim_data):
		"""
		Determine which genes are regulated by ppGpp and store the fold
		change in expression associated with each RNA.

		Attributes set:
			ppgpp_regulated_genes (ndarray[str]): RNA ID of regulated
				genes (no compartment tag)
			ppgpp_fold_changes (ndarray[float]): log2 fold change for
				each gene in ppgpp_regulated_genes
			_ppgpp_growth_parameters: parameters for interpolate.splev
				to estimate growth rate from ppGpp concentration
		"""

		def read_value(d, k):
			"""Handle empty values from raw_data as 0"""
			val = d[k]
			return 0 if val == '' else val

		# Flag to check when ppGpp regulated expression has been set
		# see set_ppgpp_expression()
		self._ppgpp_expression_set = False

		# Determine KM for ppGpp binding to RNAP
		self._solve_ppgpp_km(raw_data, sim_data)

		# Read regulation data from raw_data
		# Treats ppGpp and DksA-ppGpp regulation the same
		gene_to_rna = {g['symbol']: g['rna_id'] for g in raw_data.genes}
		regulation = {}
		for reg in raw_data.ppgpp_regulation:
			# Convert to regulated RNA
			gene = reg['Gene']
			rna = gene_to_rna.get(gene, None)
			if rna is None:
				continue

			# Add additional gene symbols for matching FC data
			curated_gene = reg['Curated Gene']
			if gene != curated_gene:
				gene_to_rna[curated_gene] = rna

			# Update value (some genes are repeated in raw_data))
			direction = read_value(reg, 'ppGpp') + read_value(reg, 'DksA-ppGpp')
			regulation[rna] = regulation.get(rna, 0) + direction

		# Read fold change data from raw_data
		## Categories A-D are statistically significant fold changes
		## Categories E-G are not significant or data is not usable
		valid_categories = {'A', 'B', 'C', 'D'}
		sample_time = 5  # Could also be 10 (5 min minimizes downstream regulation impacts)
		sample_id = '1+2+ {} min'.format(sample_time)  # Column contains FC data for given time
		rna_fold_changes = {}
		for fc in raw_data.ppgpp_fc:
			# Convert to regulated RNA
			gene = fc['Gene']
			rna = gene_to_rna.get(gene, None)
			if rna is None:
				continue

			category = fc['{} Category'.format(sample_id)]
			if category not in valid_categories:
				continue
			rna_fold_changes[rna] = fc[sample_id]

		# Store arrays of regulation
		regulated_genes = []
		regulation_direction = []
		fold_changes = []
		for rna in sorted(regulation):
			reg_dir = regulation[rna]
			fc_dir = rna_fold_changes.get(rna, 0)

			# Ignore inconsistent regulatory directions
			if reg_dir == 0 or reg_dir * fc_dir < 0:
				continue

			regulated_genes.append(rna)
			regulation_direction.append(np.sign(reg_dir))
			fold_changes.append(fc_dir)
		self.ppgpp_regulated_genes = np.array(regulated_genes)
		regulation_direction = np.array(regulation_direction)

		# Replace fold changes without data with the average
		fold_changes = np.array(fold_changes)
		average_positive_fc = fold_changes[fold_changes > 0].mean()
		fold_changes[(fold_changes == 0) & (regulation_direction < 0)] = self._fit_ppgpp_fc
		fold_changes[(fold_changes == 0) & (regulation_direction > 0)] = average_positive_fc
		self.ppgpp_fold_changes = fold_changes

		# Predict growth rate from ppGpp level
		per_dry_mass_to_per_volume = sim_data.constants.cell_density * sim_data.mass.cell_dry_mass_fraction
		ppgpp = np.array([(d['ppGpp_conc'] * per_dry_mass_to_per_volume).asNumber(PPGPP_CONC_UNITS)
			for d in raw_data.growth_rate_dependent_parameters])
		growth_rates = np.log(2) / np.array([d['doublingTime'].asNumber(units.s)
			for d in raw_data.growth_rate_dependent_parameters])
		self._ppgpp_growth_parameters = interpolate.splrep(ppgpp[::-1], growth_rates[::-1], k=1)

		if PRINT_VALUES:
			print('Supplement value (KM): {:.1f}'
				.format(np.sqrt(self._ppgpp_km_squared)))
			print('Supplement value (FC): [{:.2f}, {:.2f}]'
				.format(fold_changes.min(), fold_changes.max()))
			print('Supplement value (FC-): {:.2f}'.format(self._fit_ppgpp_fc))
			print('Supplement value (FC+): {:.2f}'.format(average_positive_fc))

	def _build_rna_data(self, raw_data, sim_data):
		"""
		Build RNA-associated simulation data from raw data.
		"""
		self._basal_rna_fractions = sim_data.mass.get_basal_rna_fractions()

		# Filter out RNAs without sequences
		all_rnas = [
			rna for rna in raw_data.rnas
			if sim_data.getter.is_valid_molecule(rna['id'])]

		# Load RNA IDs with compartment tags
		rna_ids = [rna['id'] for rna in all_rnas]
		compartments = sim_data.getter.get_compartments(rna_ids)

		rna_ids_with_compartments = [
			f'{rna_id}[{loc[0]}]' for (rna_id, loc)
			in zip(rna_ids, compartments)]

		# Load set of mRNA ids
		mRNA_ids = set([rna['id'] for rna in all_rnas if rna['type'] == 'mRNA'])

		# Load RNA half lives
		rna_id_to_half_life = {}
		reported_mRNA_half_lives = []

		for rna in raw_data.rna_half_lives:
			rna_id_to_half_life[rna['id']] = rna['half_life']

			if rna['id'] in mRNA_ids:
				reported_mRNA_half_lives.append(rna['half_life'])

		# Calculate average reported half lives of mRNAs
		average_mRNA_half_lives = np.array(reported_mRNA_half_lives).mean()

		# Get half life of each RNA - if the half life is not given, use the
		# average reported half life of mRNAs
		half_lives = np.array([
			rna_id_to_half_life.get(rna['id'], average_mRNA_half_lives).asNumber(units.s)
			for rna in all_rnas])

		# Convert to degradation rates
		rna_deg_rates = np.log(2) / half_lives

		# Load RNA expression from RNA-seq data
		expression = []
		rna_id_to_gene_id = {
			gene['rna_id']: gene['id'] for gene in raw_data.genes}

		seq_data = {
			x['Gene']: x[sim_data.basal_expression_condition]
			for x in getattr(raw_data.rna_seq_data, f'rnaseq_{RNA_SEQ_ANALYSIS}_mean')}

		for rna in all_rnas:
			gene_id = rna_id_to_gene_id[rna['id']]
			# If sequencing data is not found, initialize expression to zero.
			expression.append(seq_data.get(gene_id, 0.))

		expression = np.array(expression)

		# Calculate synthesis probabilities from expression and normalize
		synth_prob = expression*(
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
			+ rna_deg_rates
			)
		synth_prob /= synth_prob.sum()

		# Calculate EndoRNase Km values
		Km = (KCAT_ENDO_RNASE*ESTIMATE_ENDO_RNASES/rna_deg_rates) - expression

		# Load gene IDs
		gene_ids = np.array(
			[rna_id_to_gene_id[rna['id']] for rna in all_rnas])

		# Construct boolean arrays for ribosomal protein and RNAP genes
		n_rnas = len(rna_ids_with_compartments)
		is_ribosomal_protein = np.zeros(n_rnas, dtype=np.bool)
		is_RNAP	= np.zeros(n_rnas, dtype=np.bool)

		for i, rna in enumerate(all_rnas):
			for monomer_id in rna['monomer_ids']:
				if monomer_id + '[c]' in sim_data.molecule_groups.ribosomal_proteins:
					is_ribosomal_protein[i] = True
				if monomer_id + '[c]' in sim_data.molecule_groups.RNAP_subunits:
					is_RNAP[i] = True

		# Construct boolean arrays and index arrays for each rRNA type
		is_23S = np.zeros(n_rnas, dtype = np.bool)
		is_16S = np.zeros(n_rnas, dtype = np.bool)
		is_5S = np.zeros(n_rnas, dtype = np.bool)
		idx_23S = []
		idx_16S = []
		idx_5S = []

		for rnaIndex, rna in enumerate(all_rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is_23S[rnaIndex] = True
				idx_23S.append(rnaIndex)

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is_16S[rnaIndex] = True
				idx_16S.append(rnaIndex)

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is_5S[rnaIndex] = True
				idx_5S.append(rnaIndex)

		idx_23S = np.array(idx_23S)
		idx_16S = np.array(idx_16S)
		idx_5S = np.array(idx_5S)

		# Load RNA sequences and molecular weights from getter functions
		rna_seqs = sim_data.getter.get_sequences(rna_ids)
		mws = sim_data.getter.get_masses(rna_ids).asNumber(units.g / units.mol)

		# Calculate lengths and nt counts from sequence
		rna_lengths = np.array([len(seq) for seq in rna_seqs])

		# Get RNA nucleotide compositions
		ntp_abbreviations = [ntp_id[0] for ntp_id in sim_data.molecule_groups.ntps]
		nt_counts = []
		for seq in rna_seqs:
			nt_counts.append(
				[seq.count(letter) for letter in ntp_abbreviations])
		nt_counts = np.array(nt_counts)

		# Get index of gene corresponding to each RNA
		rna_id_to_gene_index = {gene['rna_id']: i
			for i, gene in enumerate(raw_data.genes)}

		# Get list of coordinates and directions for each gene
		coordinate_list = [
			gene["left_end_pos"] if gene["direction"] == '+'
			else gene["right_end_pos"]
			for gene in raw_data.genes]
		direction_list = [gene["direction"] for gene in raw_data.genes]

		# Get coordinates of oriC and terC
		oric_coordinate = sim_data.constants.oriC_center.asNumber()
		terc_coordinate = sim_data.constants.terC_center.asNumber()
		genome_length = len(raw_data.genome_sequence)

		def get_relative_coordinates(coordinates):
			"""
			Returns the genomic coordinates of a given gene coordinate relative
			to the origin of replication.
			"""
			if coordinates < terc_coordinate:
				relative_coordinates = genome_length - oric_coordinate + coordinates
			elif coordinates < oric_coordinate:
				relative_coordinates = coordinates - oric_coordinate + 1
			else:
				relative_coordinates = coordinates - oric_coordinate

			return relative_coordinates

		# Get location of transcription initiation relative to origin
		replication_coordinate = [
			get_relative_coordinates(coordinate_list[rna_id_to_gene_index[rna["id"]]])
			for rna in all_rnas]

		# Get direction of transcription
		direction = [
			(direction_list[rna_id_to_gene_index[rna["id"]]] == "+")
			for rna in all_rnas]

		# Set the lengths, nucleotide counts, molecular weights, and sequences
		# of each type of rRNAs to be identical to those of the first rRNA
		# operon. Later in the sim, transcription of all rRNA genes are set to
		# produce the rRNAs of the first operon. This is done to simplify the
		# complexation reactions that form ribosomes. In reality, all of these
		# genes produce rRNA molecules with slightly different sequences and
		# molecular weights.
		rna_lengths[idx_23S] = rna_lengths[idx_23S[0]]
		rna_lengths[idx_16S] = rna_lengths[idx_16S[0]]
		rna_lengths[idx_5S] = rna_lengths[idx_5S[0]]

		nt_counts[idx_23S, :] = nt_counts[idx_23S[0], :]
		nt_counts[idx_16S, :] = nt_counts[idx_16S[0], :]
		nt_counts[idx_5S, :] = nt_counts[idx_5S[0], :]

		mws[idx_23S] = mws[idx_23S[0]]
		mws[idx_16S] = mws[idx_16S[0]]
		mws[idx_5S] = mws[idx_5S[0]]

		id_length = max(len(id_) for id_ in rna_ids_with_compartments)
		gene_id_length = max(len(id_) for id_ in gene_ids)
		rna_data = np.zeros(
			n_rnas,
			dtype = [
				('id', 'U{}'.format(id_length)),
				('deg_rate', 'f8'),
				('length', 'i8'),
				('counts_ACGU', '4i8'),
				('mw', 'f8'),
				('is_mRNA', 'bool'),
				('is_miscRNA', 'bool'),
				('is_rRNA', 'bool'),
				('is_tRNA', 'bool'),
				('is_23S_rRNA', 'bool'),
				('is_16S_rRNA', 'bool'),
				('is_5S_rRNA', 'bool'),
				('is_ribosomal_protein', 'bool'),
				('is_RNAP',	'bool'),
				('gene_id', 'U{}'.format(gene_id_length)),
				('Km_endoRNase', 'f8'),
				('replication_coordinate', 'int64'),
				('direction', 'bool'),
				]
			)

		rna_data['id'] = rna_ids_with_compartments
		rna_data['deg_rate'] = rna_deg_rates
		rna_data['length'] = rna_lengths
		rna_data['counts_ACGU'] = nt_counts
		rna_data['mw'] = mws
		rna_data['is_mRNA'] = [rna["type"] == "mRNA" for rna in all_rnas]
		rna_data['is_miscRNA'] = [rna["type"] == "miscRNA" for rna in all_rnas]
		rna_data['is_rRNA'] = [rna["type"] == "rRNA" for rna in all_rnas]
		rna_data['is_tRNA'] = [rna["type"] == "tRNA" for rna in all_rnas]
		rna_data['is_ribosomal_protein'] = is_ribosomal_protein
		rna_data['is_RNAP'] = is_RNAP
		rna_data['is_23S_rRNA'] = is_23S
		rna_data['is_16S_rRNA'] = is_16S
		rna_data['is_5S_rRNA'] = is_5S
		rna_data['gene_id'] = gene_ids
		rna_data['Km_endoRNase'] = Km
		rna_data['replication_coordinate'] = replication_coordinate
		rna_data['direction'] = direction

		field_units = {
			'id': None,
			'deg_rate': 1 / units.s,
			'length': units.nt,
			'counts_ACGU': units.nt,
			'mw': units.g / units.mol,
			'is_mRNA': None,
			'is_miscRNA': None,
			'is_rRNA': None,
			'is_tRNA': None,
			'is_23S_rRNA': None,
			'is_16S_rRNA': None,
			'is_5S_rRNA': None,
			'is_ribosomal_protein': None,
			'is_RNAP': None,
			'gene_id': None,
			'Km_endoRNase': units.mol / units.L,
			'replication_coordinate': None,
			'direction': None,
			}

		self.rna_expression = {}
		self.rna_synth_prob = {}

		# Set basal expression and synthesis probabilities - conditional values
        # are set in the parca.
		self.rna_expression["basal"] = expression / expression.sum()
		self.rna_synth_prob["basal"] = synth_prob / synth_prob.sum()

		self.rna_data = UnitStructArray(rna_data, field_units)

	def _build_transcription(self, raw_data, sim_data):
		"""
		Build transcription-associated simulation data from raw data.
		"""
		# Load sequence data
		rna_seqs = sim_data.getter.get_sequences(
			[rna_id[:-3] for rna_id in self.rna_data['id']])

		rrna_types = ['is_23S_rRNA', 'is_16S_rRNA', 'is_5S_rRNA']
		for rrna in rrna_types:
			rrna_idx = np.where(self.rna_data[rrna])[0]
			for idx in rrna_idx[1:]:
				rna_seqs[idx] = rna_seqs[rrna_idx[0]]

		# Construct transcription sequence matrix
		maxLen = np.int64(
			self.rna_data["length"].asNumber().max()
			+ self.max_time_step * sim_data.growth_rate_parameters.rnaPolymeraseElongationRate.asNumber(units.nt / units.s)
			)

		self.transcription_sequences = np.full((len(rna_seqs), maxLen), polymerize.PAD_VALUE, dtype=np.int8)
		ntMapping = {ntpId: i for i, ntpId in enumerate(["A", "C", "G", "U"])}
		for i, sequence in enumerate(rna_seqs):
			for j, letter in enumerate(sequence):
				self.transcription_sequences[i, j] = ntMapping[letter]

		# Calculate weights of transcript nucleotide monomers
		self.transcription_monomer_weights = (
			(
				sim_data.getter.get_masses(sim_data.molecule_groups.ntps)
				- sim_data.getter.get_masses([sim_data.molecule_ids.ppi])
				)
			/ sim_data.constants.n_avogadro
			).asNumber(units.fg)

		self.transcription_end_weight = ((sim_data.getter.get_masses([sim_data.molecule_ids.ppi])
			/ sim_data.constants.n_avogadro).asNumber(units.fg))

	def _build_charged_trna(self, raw_data, sim_data):
		'''
		Loads information and creates data structures necessary for charging of tRNA

		Note:
			Requires self.rnaData so can't be built in translation even if some
			data structures would be more appropriate there.
		'''

		# Create list of charged tRNAs
		trna_names = self.rna_data['id'][self.rna_data['is_tRNA']]
		charged_trnas = [x['modified_forms'] for x in raw_data.rnas if x['id'] + '[c]' in trna_names]
		filtered_charged_trna = []
		for charged_list in charged_trnas:
			for trna in charged_list:
				# Skip modified forms so only one charged tRNA per uncharged tRNA
				if 'FMET' in trna or 'modified' in trna:
					continue

				assert('c' in sim_data.getter.get_compartment(trna))
				filtered_charged_trna += [trna + '[c]']

		self.charged_trna_names = filtered_charged_trna
		assert(len(self.charged_trna_names) == len(trna_names))

		# Create mapping of each tRNA/charged tRNA to associated AA
		trna_dict = {
			'RNA0-300[c]': 'VAL',
			'RNA0-301[c]': 'LYS',
			'RNA0-302[c]': 'LYS',
			'RNA0-303[c]': 'LYS',
			'RNA0-304[c]': 'ASN',
			'RNA0-305[c]': 'ILE',
			'RNA0-306[c]': 'MET',
			}
		aa_names = sim_data.molecule_groups.amino_acids
		aa_indices = {aa: i for i, aa in enumerate(aa_names)}
		trna_indices = {trna: i for i, trna in enumerate(trna_names)}
		self.aa_from_trna = np.zeros((len(aa_names), len(trna_names)))
		for trna in trna_names:
			aa = trna[:3].upper()
			if aa == 'ALA':
				aa = 'L-ALPHA-ALANINE'
			elif aa == 'ASP':
				aa = 'L-ASPARTATE'
			elif aa == 'SEL':
				aa = 'L-SELENOCYSTEINE'
			elif aa == 'RNA':
				aa = trna_dict[trna]

			assert('c' in sim_data.getter.get_compartment(aa))
			aa += '[c]'
			if aa in aa_names:
				aa_idx = aa_indices[aa]
				trna_idx = trna_indices[trna]
				self.aa_from_trna[aa_idx, trna_idx] = 1

		# Arrays for stoichiometry and synthetase mapping matrices
		molecules = []

		# Sparse matrix representation - i, j are row/column indices and v is value
		stoich_matrix_i = []
		stoich_matrix_j = []
		stoich_matrix_v = []

		synthetase_names = []
		synthetase_mapping_aa = []
		synthetase_mapping_syn = []

		# Get IDs of charging reactions that should be removed
		removed_reaction_ids = {
			rxn['id'] for rxn in raw_data.trna_charging_reactions_removed}

		# Get IDs of all metabolites
		metabolite_ids = {met['id'] for met in raw_data.metabolites}

		# Create stoichiometry matrix for charging reactions
		for reaction in raw_data.trna_charging_reactions:
			if reaction['id'] in removed_reaction_ids:
				continue

			# Get uncharged tRNA name for the given reaction
			trna = None
			for mol_id in reaction['stoichiometry'].keys():
				if f'{mol_id}[c]' in trna_names:
					trna = f'{mol_id}[c]'
					break

			if trna is None:
				continue

			trna_index = trna_indices[trna]

			# Get molecule information
			aa_idx = None
			for mol_id, coeff in reaction['stoichiometry'].items():

				if mol_id in metabolite_ids:
					molecule_name = "{}[{}]".format(
						mol_id, 'c'
						# Assume all metabolites are in cytosol
						)
				else:
					molecule_name = "{}[{}]".format(
						mol_id,
						sim_data.getter.get_compartment(mol_id)[0]
						)

				if molecule_name not in molecules:
					molecules.append(molecule_name)
					molecule_index = len(molecules) - 1
				else:
					molecule_index = molecules.index(molecule_name)

				aa_idx = aa_indices.get(molecule_name, aa_idx)

				# Assume coefficents given as null are -1
				if coeff is None:
					coeff = -1
				assert coeff % 1 == 0

				stoich_matrix_i.append(molecule_index)
				stoich_matrix_j.append(trna_index)
				stoich_matrix_v.append(coeff)

			assert aa_idx is not None

			# Create mapping for synthetases catalyzing charging
			for synthetase in reaction['catalyzed_by']:
				synthetase = '{}[{}]'.format(
					synthetase, sim_data.getter.get_compartment(synthetase)[0])

				if synthetase not in synthetase_names:
					synthetase_names.append(synthetase)

				synthetase_mapping_aa.append(aa_idx)
				synthetase_mapping_syn.append(synthetase_names.index(synthetase))

		# Save matrices and related lists of names
		self._stoich_matrix_i = np.array(stoich_matrix_i)
		self._stoich_matrix_j = np.array(stoich_matrix_j)
		self._stoich_matrix_v = np.array(stoich_matrix_v)

		self.aa_from_synthetase = np.zeros((len(aa_names), len(synthetase_names)))
		self.aa_from_synthetase[synthetase_mapping_aa, synthetase_mapping_syn] = 1

		self.synthetase_names = synthetase_names
		self.charging_molecules = molecules

	def charging_stoich_matrix(self):
		'''
		Creates stoich matrix from i, j, v arrays

		Returns 2D array with rows of metabolites for each tRNA charging reaction on the column
		'''

		shape = (self._stoich_matrix_i.max() + 1, self._stoich_matrix_j.max() + 1)

		out = np.zeros(shape, np.float64)
		out[self._stoich_matrix_i, self._stoich_matrix_j] = self._stoich_matrix_v

		return out

	def _build_attenuation(self, raw_data, sim_data):
		"""
		Load fold changes related to transcriptional attenuation.
		"""

		# Load data from file
		aa_trnas = []
		attenuated_rnas = []
		fold_changes = []
		gene_to_rna = {g['symbol']: g['rna_id'] for g in raw_data.genes}
		for row in raw_data.transcriptional_attenuation:
			trna_aa = row['tRNA'].split('-')[1].upper() + '[c]'
			gene = row['Target']
			rna = gene_to_rna[gene] + '[c]'

			aa_trnas.append(trna_aa)
			attenuated_rnas.append(rna)
			fold_changes.append(2**row['log2 FC'])

		self.attenuated_rna_ids = np.unique(attenuated_rnas)

		# Convert data to matrix mapping tRNA to genes with a fold change
		trna_to_row = {t: i for i, t in enumerate(sim_data.molecule_groups.amino_acids)}
		rna_to_col = {r: i for i, r in enumerate(self.attenuated_rna_ids)}
		n_aas = len(sim_data.molecule_groups.amino_acids)
		n_rnas = len(self.attenuated_rna_ids)
		self._attenuation_fold_changes = np.ones((n_aas, n_rnas))
		for trna, rna, fc in zip(aa_trnas, attenuated_rnas, fold_changes):
			i = trna_to_row[trna]
			j = rna_to_col[rna]
			self._attenuation_fold_changes[i, j] = fc

		# Attenuated RNA index mapping
		rna_to_index = {r: i for i, r in enumerate(self.rna_data['id'])}
		self.attenuated_rna_indices = np.array([rna_to_index[r] for r in self.attenuated_rna_ids])

		# Specify location in gene where attenuation will occur
		# Currently just assumes before a transcript begins elongation (position < 1)
		# TODO: base this on specific locations for each gene
		locations = np.ones(len(self.attenuated_rna_indices))
		self.attenuation_location = {idx: loc for idx, loc in zip(self.attenuated_rna_indices, locations)}

	def calculate_attenuation(self, sim_data, cell_specs):
		"""
		Calculate constants for each attenuated gene.

		TODO:
			Calculate estimate charged tRNA concentration to use instead of all tRNA
		"""

		def get_trna_conc(condition):
			spec = cell_specs[condition]
			uncharged_trna_ids = self.rna_data['id'][self.rna_data['is_tRNA']]
			counts = spec['bulkAverageContainer'].counts(uncharged_trna_ids)
			volume = (spec['avgCellDryMassInit'] / sim_data.constants.cell_density
				/ sim_data.mass.cell_dry_mass_fraction)
			# Order of operations for conc (counts last) is to get units to work well
			conc = (1 / sim_data.constants.n_avogadro / volume * counts)
			return conc

		k_units = units.umol / units.L
		trna_conc = self.aa_from_trna @ get_trna_conc('with_aa').asNumber(k_units)

		# Calculate constant for stop probability
		self.attenuation_k = np.zeros_like(self._attenuation_fold_changes)
		for i, j in zip(*np.where(self._attenuation_fold_changes != 1)):
			k = trna_conc[i] / np.log(self._attenuation_fold_changes[i, j])
			self.attenuation_k[i, j] = 1/k
		self.attenuation_k = 1 / k_units * self.attenuation_k

		# Adjust basal synthesis probabilities to account for less synthesis
		# due to attenuation
		condition = 'basal'
		basal_prob = sim_data.process.transcription_regulation.basal_prob
		delta_prob = sim_data.process.transcription_regulation.get_delta_prob_matrix()
		p_promoter_bound = np.array([
			sim_data.pPromoterBound[condition][tf]
			for tf in sim_data.process.transcription_regulation.tf_ids
			])
		delta = delta_prob @ p_promoter_bound
		basal_stop_prob = self.get_attenuation_stop_probabilities(get_trna_conc(condition))
		basal_synth_prob = (basal_prob + delta)[self.attenuated_rna_indices]
		self.attenuation_basal_prob_adjustments = basal_synth_prob * (1 / (1 - basal_stop_prob) - 1)

	def get_attenuation_stop_probabilities(self, trna_conc):
		trna_by_aa = units.matmul(self.aa_from_trna, trna_conc)
		stop_prob = 1 - np.exp(units.strip_empty_units(trna_by_aa @ self.attenuation_k))
		return stop_prob

	def _build_elongation_rates(self, raw_data, sim_data):
		self.max_elongation_rate = sim_data.constants.RNAP_elongation_rate_max.asNumber(units.nt / units.s)
		self.rRNA_indexes = np.where(self.rna_data['is_rRNA'])[0]

	def make_elongation_rates(self, random, base, time_step, variable_elongation=False):
		return make_elongation_rates(
			random,
			self.transcription_sequences.shape[0],
			base,
			self.rRNA_indexes,
			self.max_elongation_rate,
			time_step,
			variable_elongation)

	def _solve_ppgpp_km(self, raw_data, sim_data):
		"""
		Solves for general expression rates for bound and free RNAP and
		a KM for ppGpp to RNAP based on global cellular measurements.
		Parameters are solved for at different doubling times using a
		gradient descent method to minimize the difference in expression
		of stable RNA compared to the measured RNA in a cell.
		Assumes a Hill coefficient of 2 for ppGpp binding to RNAP.

		Attributes set:
			_fit_ppgpp_fc (float): log2 fold change in stable RNA expression
				from a fast doubling time to a slow doubling time based on
				the rates of bound and free RNAP expression found
			_ppgpp_km_squared (float): squared and unitless KM value for
				to limit computation needed for fraction bound
			ppgpp_km (float with mol / volume units): KM for ppGpp binding
				to RNAP
		"""

		# Data for different doubling times (100, 60, 40, 30, 24 min)
		per_dry_mass_to_per_volume = sim_data.constants.cell_density * sim_data.mass.cell_dry_mass_fraction
		ppgpp = np.array(
			[(d['ppGpp_conc'] * per_dry_mass_to_per_volume).asNumber(PPGPP_CONC_UNITS)
			for d in raw_data.growth_rate_dependent_parameters])**2
		rna = np.array([d['rnaMassFraction']
			for d in raw_data.dry_mass_composition])
		mass_per_cell = np.array([d['averageDryMass'].asNumber(units.fg)
			for d in raw_data.dry_mass_composition])
		rnap_per_cell = np.array([d['RNAP_per_cell']
			for d in raw_data.growth_rate_dependent_parameters])

		# Variables for the objective
		## a1: rate of RNA production from free RNAP
		## a2: rate of RNA production from RNAP bound to ppGpp
		## km: KM of ppGpp binding to RNAP
		a1s, a2s, kms = sp.symbols('a1 a2 km')

		# Create objective to minimize
		## Objective is squared difference between RNA created with different rates for RNAP bound
		## to ppGpp and free RNAP compared to measured RNA in the cell for each measured doubling time.
		## Use sp.exp to prevent negative parameter values, also improves stability for larger step size.
		difference = rnap_per_cell / mass_per_cell * (sp.exp(a1s)*(1 - ppgpp/(sp.exp(kms) + ppgpp))
			+ sp.exp(a2s)*ppgpp/(sp.exp(kms) + ppgpp)) - rna
		J = difference.dot(difference)

		# Convert to functions for faster performance
		dJda1 = sp.lambdify((a1s, a2s, kms), J.diff(a1s))
		dJda2 = sp.lambdify((a1s, a2s, kms), J.diff(a2s))
		dJdkm = sp.lambdify((a1s, a2s, kms), J.diff(kms))
		J = sp.lambdify((a1s, a2s, kms), J)

		# Initial parameters
		a1 = np.log(0.02)
		a2 = np.log(0.01)
		km = np.log(575)
		step_size = 1.

		# Use gradient descent to find
		obj = J(a1, a2, km)
		old_obj = 100
		step = 0
		max_step = 1e5
		tol = 1e-6
		rel_tol = 1e-9
		while obj > tol and 1 - obj / old_obj > rel_tol:
			a1 -= dJda1(a1, a2, km) * step_size
			a2 -= dJda2(a1, a2, km) * step_size
			km -= dJdkm(a1, a2, km) * step_size

			old_obj = obj
			obj = J(a1, a2, km)

			step += 1
			if step > max_step:
				raise RuntimeError('Fitting ppGpp binding KM failed to converge.'
					' Check tolerances or maximum number of steps.')

		a1 = np.exp(a1)
		a2 = np.exp(a2)
		km = np.exp(km)
		f_low = ppgpp[-1] / (km + ppgpp[-1])
		fc = np.log2(a2 / (a1*(1 - f_low) + a2*f_low))

		self._fit_ppgpp_fc = fc
		self._ppgpp_km_squared = km
		self.ppgpp_km = np.sqrt(km) * PPGPP_CONC_UNITS

	def get_rna_fractions(self, ppgpp):
		"""
		Calculates expected RNA subgroup mass fractions based on ppGpp
		concentration. If ppGpp expression has not been set yet, uses default
		measured fractions.

		Args:
			ppgpp (float with or without mol / volume units): concentration of ppGpp,
				if unitless, should represent the concentration of PPGPP_CONC_UNITS

		Returns:
			dict[str, float]: mass fraction for each subgroup mass, values sum to 1
		"""

		if self._ppgpp_expression_set:
			exp = self.expression_from_ppgpp(ppgpp)
			mass = self.rna_data['mw'] * exp
			mass = (mass / units.sum(mass)).asNumber()
			fractions =  {
				'23S': mass[self.rna_data['is_23S_rRNA']].sum(),
				'16S': mass[self.rna_data['is_16S_rRNA']].sum(),
				'5S': mass[self.rna_data['is_5S_rRNA']].sum(),
				'trna': mass[self.rna_data['is_tRNA']].sum(),
				'mrna': mass[self.rna_data['is_mRNA']].sum(),
				}
		else:
			fractions = self._basal_rna_fractions

		return fractions

	def set_ppgpp_expression(self, sim_data):
		"""
		Called during the parca to determine expression of each gene for ppGpp
		bound and free RNAP.

		Attributes set:
			exp_ppgpp (ndarray[float]): expression for each gene when RNAP
				is bound to ppGpp
			exp_free (ndarray[float]): expression for each gene when RNAP
				is not bound to ppGpp
			ppgpp_km (float with units of mol / vol): KM for ppGpp binding to RNAP
			ppgpp_km_squared (float): squared and unitless version of KM for
				faster computation in other functions
		"""

		ppgpp_aa = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time['with_aa'])
		ppgpp_basal = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time['basal'])
		f_ppgpp_aa = self.fraction_rnap_bound_ppgpp(ppgpp_aa)
		f_ppgpp_basal = self.fraction_rnap_bound_ppgpp(ppgpp_basal)

		rna_idx = {r[:-3]: i for i, r in enumerate(self.rna_data['id'])}
		fcs = np.zeros(len(self.rna_data))
		for rna, fc in zip(self.ppgpp_regulated_genes, self.ppgpp_fold_changes):
			fcs[rna_idx[rna]] = fc
		exp = self.rna_expression['basal']
		self.exp_ppgpp = ((2**fcs * exp * (1 - f_ppgpp_aa) / (1 - f_ppgpp_basal))
			/ (1 - 2**fcs * (f_ppgpp_aa - f_ppgpp_basal * (1 - f_ppgpp_aa) / (1 - f_ppgpp_basal))))
		self.exp_free = (exp - self.exp_ppgpp*f_ppgpp_basal) / (1 - f_ppgpp_basal)
		self.exp_free[self.exp_free < 0] = 0  # fold change is limited by KM, can't have very high positive fold changes

		self._ppgpp_expression_set = True

	def adjust_polymerizing_ppgpp_expression(self, sim_data):
		"""
		Adjust ppGpp expression based on fit for ribosome and RNAP physiological constraints
		using least squares fit for 3 conditions with different growth rates/ppGpp.

		Modifies attributes:
			exp_ppgpp (ndarray[float]): expression for each gene when RNAP
				is bound to ppGpp, adjusted for necessary RNAP and ribosome
				expression, normalized to 1
			exp_free (ndarray[float]): expression for each gene when RNAP
				is not bound to ppGpp, adjusted for necessary RNAP and ribosome
				expression, normalized to 1

		Note:
			See docs/processes/transcription_regulation.pdf for a description
			of the math used in this section.

		TODO:
			fit for all conditions and not just those specified below?
		"""

		# Fraction RNAP bound to ppGpp in different conditions
		ppgpp_aa = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time['with_aa'])
		ppgpp_basal = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time['basal'])
		ppgpp_anaerobic = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time['no_oxygen'])
		f_ppgpp_aa = self.fraction_rnap_bound_ppgpp(ppgpp_aa)
		f_ppgpp_basal = self.fraction_rnap_bound_ppgpp(ppgpp_basal)
		f_ppgpp_anaerobic = self.fraction_rnap_bound_ppgpp(ppgpp_anaerobic)

		# Adjustments for TFs
		tf_adjustments = {}
		delta_prob = sim_data.process.transcription_regulation.get_delta_prob_matrix()
		for condition in ['with_aa', 'basal', 'no_oxygen']:
			p_promoter_bound = np.array([
				sim_data.pPromoterBound[condition][tf]
				for tf in sim_data.process.transcription_regulation.tf_ids
				])
			delta = delta_prob @ p_promoter_bound
			# TODO(jerry): Handle the NaN or Infinity values from this line then
			#  use `with np.errstate(invalid='ignore'):` to suppress its warnings.
			tf_adjustments[condition] = delta / sim_data.process.transcription.rna_synth_prob[condition]

		# Solve least squares fit for expression of each component of RNAP and ribosomes
		self._normalize_ppgpp_expression()  # Need to normalize first to get correct scale
		adjusted_mask = self.rna_data['is_RNAP'] | self.rna_data['is_ribosomal_protein'] | self.rna_data['is_rRNA']
		F = np.array([[1- f_ppgpp_aa, f_ppgpp_aa], [1 - f_ppgpp_basal, f_ppgpp_basal], [1 - f_ppgpp_anaerobic, f_ppgpp_anaerobic]])
		Flst = np.linalg.inv(F.T.dot(F)).dot(F.T)
		expression = np.array([
			self.rna_expression['with_aa'] * (1 - tf_adjustments['with_aa']),
			self.rna_expression['basal'] * (1 - tf_adjustments['basal']),
			self.rna_expression['no_oxygen'] * (1 - tf_adjustments['no_oxygen'])])
		adjusted_free, adjusted_ppgpp = Flst.dot(expression)
		self.exp_free[adjusted_mask] = adjusted_free[adjusted_mask]
		self.exp_ppgpp[adjusted_mask] = adjusted_ppgpp[adjusted_mask]
		self._normalize_ppgpp_expression()

	def adjust_ppgpp_expression_for_tfs(self, sim_data):
		"""
		Adjusts ppGpp regulated expression to get expression with and without
		ppGpp regulation to match in basal condition and taking into account
		the effect transcription factors will have.

		TODO:
			Should this not adjust polymerizing genes (adjusted_mask in
				adjust_polymerizing_ppgpp_expression) since they have already
				been adjusted for transcription factor effects?
		"""

		condition = 'basal'

		# Current (unnormalized) probabilities from ppGpp regulation
		ppgpp_conc = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time[condition])
		old_prob, factor = self.synth_prob_from_ppgpp(ppgpp_conc,
			sim_data.process.replication.get_average_copy_number)

		# Calculate the average expected effect of TFs in basal condition
		delta_prob = sim_data.process.transcription_regulation.get_delta_prob_matrix()
		p_promoter_bound = np.array([
			sim_data.pPromoterBound[condition][tf]
			for tf in sim_data.process.transcription_regulation.tf_ids
			])
		delta = delta_prob @ p_promoter_bound

		# Calculate the required probability to match expression without ppGpp
		new_prob = normalize(self.rna_expression[condition] * factor) - delta
		new_prob[new_prob < 0] = 0
		new_prob = normalize(new_prob)

		# Determine adjustments to the current ppGpp expression to scale
		# to the expected expression
		with np.errstate(invalid='ignore'):
			adjustment = new_prob / old_prob
		adjustment[~np.isfinite(adjustment)] = 1

		# Scale free and bound expression and renormalize ppGpp regulated expression
		self.exp_free *= adjustment
		self.exp_ppgpp *= adjustment
		self._normalize_ppgpp_expression()

	def _normalize_ppgpp_expression(self):
		"""
		Normalize both free and ppGpp bound expression values to 1.
		"""

		self.exp_free /= self.exp_free.sum()
		self.exp_ppgpp /= self.exp_ppgpp.sum()

	def fraction_rnap_bound_ppgpp(self, ppgpp):
		"""
		Calculates the fraction of RNAP expected to be bound to ppGpp
		at a given concentration of ppGpp.

		Args:
			ppgpp (float with or without mol / volume units): concentration of ppGpp,
				if unitless, should represent the concentration of PPGPP_CONC_UNITS

		Returns:
			float: fraction of RNAP that will be bound to ppGpp
		"""

		if units.hasUnit(ppgpp):
			ppgpp = ppgpp.asNumber(PPGPP_CONC_UNITS)

		return ppgpp**2 / (self._ppgpp_km_squared + ppgpp**2)

	def expression_from_ppgpp(self, ppgpp):
		"""
		Calculates the expression of each gene at a given concentration of ppGpp.

		Args:
			ppgpp (float with or without mol / volume units): concentration of ppGpp,
				if unitless, should represent the concentration of PPGPP_CONC_UNITS

		Returns:
			ndarray[float]: normalized expression for each gene
		"""

		f_ppgpp = self.fraction_rnap_bound_ppgpp(ppgpp)
		return normalize(self.exp_free * (1 - f_ppgpp) + self.exp_ppgpp * f_ppgpp)

	def synth_prob_from_ppgpp(self, ppgpp, copy_number):
		"""
		Calculates the synthesis probability of each gene at a given concentration
		of ppGpp.

		Args:
			ppgpp (float with mol / volume units): concentration of ppGpp
			copy_number (Callable[float, int]): function that gives the expected copy
				number given a doubling time and gene replication coordinate

		Returns
			prob (ndarray[float]): normalized synthesis probability for each gene
			factor (ndarray[float]): factor to adjust expression to probability for each gene

		Note:
			copy_number should be sim_data.process.replication.get_average_copy_number
			but saving the function handle as a class attribute prevents pickling of sim_data
			without additional handling
		"""

		ppgpp = ppgpp.asNumber(PPGPP_CONC_UNITS)
		f_ppgpp = self.fraction_rnap_bound_ppgpp(ppgpp)

		y = interpolate.splev(ppgpp, self._ppgpp_growth_parameters)
		growth = max(cast(float, y), 0.0)
		tau = np.log(2) / growth / 60
		loss = growth + self.rna_data['deg_rate'].asNumber(1 / units.s)
		n_avg_copy = copy_number(tau, self.rna_data['replication_coordinate'])

		# Return values
		factor = loss / n_avg_copy
		prob = normalize((self.exp_free * (1 - f_ppgpp) + self.exp_ppgpp * f_ppgpp) * factor)

		return prob, factor
