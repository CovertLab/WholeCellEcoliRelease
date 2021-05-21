"""
KnowledgeBase for Ecoli
Whole-cell knowledge base for Ecoli. Contains all raw, un-fit data processed
directly from CSV flat files.

"""
from __future__ import absolute_import, division, print_function

import io
import os
import json

from reconstruction.spreadsheets import read_tsv
from wholecell.io import tsv
from wholecell.utils import units  # used by eval()


FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = (
	"amino_acid_pathways.tsv",
	"biomass.tsv",
	"compartments.tsv",
	"complexation_reactions.tsv",
	"complexation_reactions_removed.tsv",
	"disabled_kinetic_reactions.tsv",
	"dry_mass_composition.tsv",
	"endoRNases.tsv",
	"equilibrium_reaction_rates.tsv",
	"equilibrium_reactions.tsv",
	"equilibrium_reactions_removed.tsv",
	"fold_changes.tsv",
	"fold_changes_nca.tsv",
	"fold_changes_removed.tsv",
	"genes.tsv",
	"growth_rate_dependent_parameters.tsv",
	"linked_metabolites.tsv",
	"metabolic_reactions.tsv",
	"metabolic_reactions_removed.tsv",
	"metabolism_kinetics.tsv",
	"metabolite_concentrations.tsv",
	"metabolites.tsv",
	"modified_proteins.tsv",
	"molecular_weight_keys.tsv",
	"ppgpp_fc.tsv",
	"ppgpp_regulation.tsv",
	"protein_half_lives_measured.tsv",
	"protein_half_lives_n_end_rule.tsv",
	"protein_modification_reactions.tsv",
	"proteins.tsv",
	"relative_metabolite_concentrations.tsv",
	"rna_half_lives.tsv",
	"rnas.tsv",
	"secretions.tsv",
	"sequence_motifs.tsv",
	"transcription_factors.tsv",
	"transcriptional_attenuation.tsv",
	"tf_one_component_bound.tsv",
	"translation_efficiency.tsv",
	"trna_charging_reactions.tsv",
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
	os.path.join("common_names", "genes.tsv"),
	os.path.join("common_names", "metabolites.tsv"),
	os.path.join("common_names", "proteins.tsv"),
	os.path.join("common_names", "reactions.tsv"),
	os.path.join("common_names", "rnas.tsv"),
	os.path.join("base_codes", "amino_acids.tsv"),
	os.path.join("base_codes", "ntp.tsv"),
	os.path.join("base_codes", "dntp.tsv"),
	os.path.join("adjustments", "translation_efficiencies_adjustments.tsv"),
	os.path.join("adjustments", "rna_expression_adjustments.tsv"),
	os.path.join("adjustments", "rna_deg_rates_adjustments.tsv"),
	os.path.join("adjustments", "protein_deg_rates_adjustments.tsv"),
	os.path.join("adjustments", "relative_metabolite_concentrations_changes.tsv"),
	)
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = (
	"parameters.tsv",
	"mass_parameters.tsv",
	"dna_supercoiling.tsv"
	)

class DataStore(object):
	def __init__(self):
		pass

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self):
		# Load raw data from TSV files
		for filename in LIST_OF_DICT_FILENAMES:
			self._load_tsv(FLAT_DIR, os.path.join(FLAT_DIR, filename))

		for filename in LIST_OF_PARAMETER_FILENAMES:
			self._load_parameters(os.path.join(FLAT_DIR, filename))

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
