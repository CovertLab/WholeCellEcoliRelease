"""
KnowledgeBase for Ecoli
Whole-cell knowledge base for Ecoli. Contains all raw, un-fit data processed
directly from CSV flat files.

"""
from __future__ import absolute_import, division, print_function

import io
import os
import json
from typing import List

from reconstruction.spreadsheets import read_tsv
from wholecell.io import tsv
from wholecell.utils import units  # used by eval()

FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = [
	"amino_acid_export_kms.tsv",
	"amino_acid_pathways.tsv",
	"biomass.tsv",
	"compartments.tsv",
	"complexation_reactions.tsv",
	"complexation_reactions_added.tsv",
	"complexation_reactions_modified.tsv",
	"complexation_reactions_removed.tsv",
	"disabled_kinetic_reactions.tsv",
	"dna_sites.tsv",
	"dry_mass_composition.tsv",
	"endoRNases.tsv",
	"equilibrium_reaction_rates.tsv",
	"equilibrium_reactions.tsv",
	"equilibrium_reactions_added.tsv",
	"equilibrium_reactions_removed.tsv",
	"fold_changes.tsv",
	"fold_changes_nca.tsv",
	"fold_changes_removed.tsv",
	"genes.tsv",
	"growth_rate_dependent_parameters.tsv",
	"linked_metabolites.tsv",
	"metabolic_reactions.tsv",
	"metabolic_reactions_added.tsv",
	"metabolic_reactions_modified.tsv",
	"metabolic_reactions_removed.tsv",
	"metabolism_kinetics.tsv",
	"metabolite_concentrations.tsv",
	"metabolite_concentrations_removed.tsv",
	"metabolites.tsv",
	"metabolites_added.tsv",
	"modified_proteins.tsv",
	"molecular_weight_keys.tsv",
	"ppgpp_fc.tsv",
	"ppgpp_regulation.tsv",
	"ppgpp_regulation_removed.tsv",
	"protein_half_lives_measured.tsv",
	"protein_half_lives_n_end_rule.tsv",
	"proteins.tsv",
	"relative_metabolite_concentrations.tsv",
	"rna_half_lives.tsv",
	"rnas.tsv",
	"secretions.tsv",
	"sequence_motifs.tsv",
	"transcription_factors.tsv",
	# "transcription_units.tsv",  # special cased in the constructor
	"transcription_units_modified.tsv",
	"transcription_units_removed.tsv",
	"transcriptional_attenuation.tsv",
	"transcriptional_attenuation_removed.tsv",
	"tf_one_component_bound.tsv",
	"translation_efficiency.tsv",
	"trna_charging_reactions.tsv",
	"trna_charging_reactions_added.tsv",
	"trna_charging_reactions_removed.tsv",
	"two_component_systems.tsv",
	"two_component_system_templates.tsv",
	os.path.join("mass_fractions", "glycogen_fractions.tsv"),
	os.path.join("mass_fractions", "ion_fractions.tsv"),
	os.path.join("mass_fractions", "LPS_fractions.tsv"),
	os.path.join("mass_fractions", "lipid_fractions.tsv"),
	os.path.join("mass_fractions", "murein_fractions.tsv"),
	os.path.join("mass_fractions", "soluble_fractions.tsv"),
	os.path.join("trna_data","trna_ratio_to_16SrRNA_0p4.tsv"),
	os.path.join("trna_data","trna_ratio_to_16SrRNA_0p7.tsv"),
	os.path.join("trna_data","trna_ratio_to_16SrRNA_1p6.tsv"),
	os.path.join("trna_data","trna_ratio_to_16SrRNA_1p07.tsv"),
	os.path.join("trna_data","trna_ratio_to_16SrRNA_2p5.tsv"),
	os.path.join("trna_data","trna_growth_rates.tsv"),
	os.path.join("rna_seq_data","rnaseq_rsem_tpm_mean.tsv"),
	os.path.join("rna_seq_data","rnaseq_rsem_tpm_std.tsv"),
	os.path.join("rna_seq_data","rnaseq_seal_rpkm_mean.tsv"),
	os.path.join("rna_seq_data","rnaseq_seal_rpkm_std.tsv"),
	os.path.join("condition", "tf_condition.tsv"),
	os.path.join("condition", "condition_defs.tsv"),
	os.path.join("condition", "environment_molecules.tsv"),
	os.path.join("condition", "timelines_def.tsv"),
	os.path.join("condition", "media_recipes.tsv"),
	os.path.join("condition", "media", "M9.tsv"),
	os.path.join("condition", "media", "M9_GLC.tsv"),
	os.path.join("condition", "media", "5X_supplement_EZ.tsv"),
	os.path.join("base_codes", "amino_acids.tsv"),
	os.path.join("base_codes", "ntp.tsv"),
	os.path.join("base_codes", "dntp.tsv"),
	os.path.join("adjustments", "amino_acid_pathways.tsv"),
	os.path.join("adjustments", "translation_efficiencies_adjustments.tsv"),
	os.path.join("adjustments", "rna_expression_adjustments.tsv"),
	os.path.join("adjustments", "rna_deg_rates_adjustments.tsv"),
	os.path.join("adjustments", "protein_deg_rates_adjustments.tsv"),
	os.path.join("adjustments", "relative_metabolite_concentrations_changes.tsv"),
	]
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = (
	"dna_supercoiling.tsv",
	"parameters.tsv",
	"mass_parameters.tsv",
	)

REMOVED_DATA = {
	'complexation_reactions': 'complexation_reactions_removed',
	'equilibrium_reactions': 'equilibrium_reactions_removed',
	'fold_changes': 'fold_changes_removed',
	'fold_changes_nca': 'fold_changes_removed',
	'metabolic_reactions': 'metabolic_reactions_removed',
	'metabolite_concentrations': 'metabolite_concentrations_removed',
	'ppgpp_regulation': 'ppgpp_regulation_removed',
	'transcriptional_attenuation': 'transcriptional_attenuation_removed',
	'trna_charging_reactions': 'trna_charging_reactions_removed',
	}
MODIFIED_DATA = {
	'complexation_reactions': 'complexation_reactions_modified',
	'metabolic_reactions': 'metabolic_reactions_modified',
	}

ADDED_DATA = {
	'complexation_reactions': 'complexation_reactions_added',
	'equilibrium_reactions': 'equilibrium_reactions_added',
	'metabolic_reactions': 'metabolic_reactions_added',
	'metabolites': 'metabolites_added',
	'trna_charging_reactions': 'trna_charging_reactions_added',
	}

class DataStore(object):
	def __init__(self):
		pass

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self, operons_on: bool):
		self.operons_on = operons_on

		self.compartments: List[dict] = []  # mypy can't track setattr(self, attr_name, rows)
		self.transcription_units: List[dict] = []

		if self.operons_on:
			LIST_OF_DICT_FILENAMES.append('transcription_units.tsv')
			REMOVED_DATA.update({
				'transcription_units': 'transcription_units_removed',
				})
			MODIFIED_DATA.update({
				'transcription_units': 'transcription_units_modified',
				})

		# Load raw data from TSV files
		for filename in LIST_OF_DICT_FILENAMES:
			self._load_tsv(FLAT_DIR, os.path.join(FLAT_DIR, filename))

		for filename in LIST_OF_PARAMETER_FILENAMES:
			self._load_parameters(os.path.join(FLAT_DIR, filename))

		self._prune_data()
		self._join_data()
		self._modify_data()

		self.genome_sequence = self._load_sequence(os.path.join(FLAT_DIR, SEQUENCE_FILE))

	def _load_tsv(self, dir_name, file_name):
		path = self
		for sub_path in file_name[len(dir_name) + 1 : ].split(os.path.sep)[:-1]:
			if not hasattr(path, sub_path):
				setattr(path, sub_path, DataStore())
			path = getattr(path, sub_path)
		attr_name = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(path, attr_name, [])

		rows = read_tsv(file_name)
		setattr(path, attr_name, rows)

	def _load_sequence(self, file_path):
		from Bio import SeqIO

		with open(file_path, "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				return record.seq

	def _load_parameters(self, file_path):
		attr_name = file_path.split(os.path.sep)[-1].split(".")[0]
		param_dict = {}

		with io.open(file_path, "rb") as csvfile:
			reader = tsv.dict_reader(csvfile)

			for row in reader:
				value = json.loads(row['value'])
				if row['units'] != '':
					# `eval()` the units [risky!] then strip it to just a unit
					# since `a_list * a_float` (like `1.0 [1/s]`) fails, and
					# `a_list * an_int` repeats the list, which is also broken.
					unit = eval(row['units'])   # risky!
					unit = units.getUnit(unit)  # strip
					value = value * unit
				param_dict[row['name']] = value

		setattr(self, attr_name, param_dict)

	def _prune_data(self):
		"""
		Remove rows that are specified to be removed. Data will only be removed
		if all data in a row in the file specifying rows to be removed matches
		the same columns in the raw data file.
		"""

		# Check each pair of files to be removed
		for data_attr, attr_to_remove in REMOVED_DATA.items():
			# Build the set of data to identify rows to be removed
			data_to_remove = getattr(self, attr_to_remove)
			removed_cols = list(data_to_remove[0].keys())
			removed_ids = set()
			for row in data_to_remove:
				removed_ids.add(tuple([row[col] for col in removed_cols]))

			# Remove any matching rows
			data = getattr(self, data_attr)
			n_entries = len(data)
			for i, row in enumerate(data[::-1]):
				checked_id = tuple([row[col] for col in removed_cols])
				if checked_id in removed_ids:
					data.pop(n_entries - i - 1)

	def _join_data(self):
		"""
		Add rows that are specified in additional files. Data will only be added
		if all the loaded columns from both datasets match.
		"""

		# Join data for each file with data to be added
		for data_attr, attr_to_add in ADDED_DATA.items():
			# Get datasets to join
			data = getattr(self, data_attr)
			added_data = getattr(self, attr_to_add)

			# Check columns are the same for each dataset
			col_diff = set(data[0].keys()).symmetric_difference(added_data[0].keys())
			if col_diff:
				raise ValueError(f'Could not join datasets {data_attr} and {attr_to_add} '
					f'because columns do not match (different columns: {col_diff}).')

			# Join datasets
			for row in added_data:
				data.append(row)

	def _modify_data(self):
		"""
		Modify entires in rows that are specified to be modified. Rows must be
		identified by their entries in the first column (usually the ID column).
		"""
		# Check each pair of files to be modified
		for data_attr, modify_attr in MODIFIED_DATA.items():
			# Build the set of data to identify rows to be modified
			data_to_modify = getattr(self, modify_attr)
			id_col_name = list(data_to_modify[0].keys())[0]

			id_to_modified_cols = {}
			for row in data_to_modify:
				id_to_modified_cols[row[id_col_name]] = row

			# Modify any matching rows with identical IDs
			data = getattr(self, data_attr)

			if list(data[0].keys())[0] != id_col_name:
				raise ValueError(f'Could not modify data {data_attr} with '
					f'{modify_attr} because the names of the first columns '
					f'do not match.')

			modified_entry_ids = set()
			for i, row in enumerate(data):
				if row[id_col_name] in id_to_modified_cols:
					data[i] = id_to_modified_cols[row[id_col_name]]
					modified_entry_ids.add(row[id_col_name])

			# Check for entries in modification data that do not exist in
			# original data
			id_diff = set(id_to_modified_cols.keys()).symmetric_difference(
				modified_entry_ids)
			if id_diff:
				raise ValueError(f'Could not modify data {data_attr} with '
					f'{modify_attr} because of one or more entries in '
					f'{modify_attr} that do not exist in {data_attr} '
					f'(nonexistent entries: {id_diff}).')
