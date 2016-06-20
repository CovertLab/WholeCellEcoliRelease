"""
KnowledgeBase for Ecoli

Whole-cell knowledge base for Ecoli. Contains all raw, un-fit data processed
directly from CSV flat files.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/11/2015
"""
from __future__ import division

import os
import csv
from reconstruction.spreadsheets import JsonReader
import json
from itertools import ifilter

from wholecell.utils import units

CSV_DIALECT = csv.excel_tab
FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = (
	"compartments.tsv",
	"complexationReactions.tsv",
	"enzymeKinetics.tsv",
	"genes.tsv",
	"metabolites.tsv",
	"metaboliteConcentrations.tsv",
	"modificationReactions.tsv",
	"modifiedRnas.tsv",
	"polymerized.tsv",
	"promoters.tsv",
	"proteinComplexes.tsv",
	"proteins.tsv",
	"reactions.tsv",
	"rnas.tsv",
	"terminators.tsv",
	"transcriptionUnits.tsv",
	"dryMassComposition.tsv",
	"biomass.tsv",
	"secretions.tsv",
	"water.tsv",
	"chromosome.tsv",
	"massAtReplicationInitiation.tsv",
	"equilibriumReactions.tsv",
	"foldChanges.tsv",
	"tfIds.tsv",
	"endoRnases.tsv",
	"translationEfficiency.tsv",
	os.path.join("massFractions", "glycogenFractions.tsv"),
	os.path.join("massFractions", "ionFractions.tsv"),
	os.path.join("massFractions", "LPSFractions.tsv"),
	os.path.join("massFractions", "lipidFractions.tsv"),
	os.path.join("massFractions", "mureinFractions.tsv"),
	os.path.join("massFractions", "solubleFractions.tsv"),
	os.path.join("trnaData","trna_ratio_to_16SrRNA_0p4.tsv"),
	os.path.join("trnaData","trna_ratio_to_16SrRNA_0p7.tsv"),
	os.path.join("trnaData","trna_ratio_to_16SrRNA_1p6.tsv"),
	os.path.join("trnaData","trna_ratio_to_16SrRNA_1p07.tsv"),
	os.path.join("trnaData","trna_ratio_to_16SrRNA_2p5.tsv"),
	os.path.join("trnaData","trna_growth_rates.tsv"),
	"growthRateDependentParameters.tsv",
	os.path.join("rna_seq_data","rnaseq_rsem_tpm_mean.tsv"),
	os.path.join("rna_seq_data","rnaseq_rsem_tpm_std.tsv"),
	os.path.join("rna_seq_data","rnaseq_seal_rpkm_mean.tsv"),
	os.path.join("rna_seq_data","rnaseq_seal_rpkm_std.tsv"),
	os.path.join("environment", "condition_doubling_time.tsv"),
	os.path.join("environment", "tf_condition.tsv"),
	os.path.join("environment", "condition_defs.tsv"),
	os.path.join("environment", "000000_wildtype", "nutrients_000000.tsv"),
	os.path.join("environment", "000001_cut_glucose", "nutrients_000000.tsv"),
	os.path.join("environment", "000001_cut_glucose", "nutrients_001200.tsv"),
	os.path.join("environment", "000002_add_aa", "nutrients_000000.tsv"),
	os.path.join("environment", "000002_add_aa", "nutrients_001200.tsv"),
	os.path.join("environment", "000003_aa", "nutrients_000000.tsv"),
	os.path.join("environment", "000004_oxygen_absent", "nutrients_000000.tsv"),
	os.path.join("environment", "000005_indole_present", "nutrients_000000.tsv"),
	os.path.join("environment", "000006_tungstate_present", "nutrients_000000.tsv"),
	os.path.join("environment", "000007_quercetin_present", "nutrients_000000.tsv"),
	os.path.join("environment", "000008_gallate_present", "nutrients_000000.tsv"),
	os.path.join("environment", "000009_succinate_carbon_source", "nutrients_000000.tsv"),
	os.path.join("environment", "000010_acetate_carbon_source", "nutrients_000000.tsv"),
	os.path.join("environment", "000011_fumarate_carbon_source", "nutrients_000000.tsv"),
	os.path.join("environment", "000012_malate_carbon_source", "nutrients_000000.tsv"),
	os.path.join("environment", "000013_nitrate_present", "nutrients_000000.tsv"),
	os.path.join("environment", "000014_nitrite_present", "nutrients_000000.tsv"),
	os.path.join("environment", "000015_calcium_absent", "nutrients_000000.tsv"),
	os.path.join("environment", "000016_magnesium_absent", "nutrients_000000.tsv"),
	os.path.join("environment", "000017_phosphate_absent", "nutrients_000000.tsv"),
	)
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = ("parameters.tsv", "mass_parameters.tsv")
CONSTANTS_FILENAME = "constants.tsv"

class DataStore(object):
	def __init__(self):
		pass

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self):
		# Load raw data from TSV files
		for filename in LIST_OF_DICT_FILENAMES:
			self._load_tsv(os.path.join(FLAT_DIR, filename))

		for filename in LIST_OF_PARAMETER_FILENAMES:
			self._load_parameters(os.path.join(FLAT_DIR, filename))
		self._load_parameters(os.path.join(FLAT_DIR, CONSTANTS_FILENAME))

		self.genome_sequence = self._load_sequence(os.path.join(FLAT_DIR, SEQUENCE_FILE))

	def _load_tsv(self, file_name):
		path = self
		for subPath in file_name[len(FLAT_DIR) + 1 : ].split(os.path.sep)[:-1]:
			if not hasattr(path, subPath):
				setattr(path, subPath, DataStore())
			path = getattr(path, subPath)
		attrName = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(path, attrName, [])

		with open(file_name, 'rU') as csvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", csvfile), # Strip comments
				dialect = CSV_DIALECT)
			setattr(path, attrName, [row for row in reader])

	def _load_sequence(self, file_path):
		from Bio import SeqIO
		with open(file_path, "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				return record.seq

	def _load_parameters(self, file_path):
		attrName = file_path.split(os.path.sep)[-1].split(".")[0]
		paramDict = {}
		with open(file_path, "rU") as csvfile:
			reader = csv.DictReader(csvfile, dialect = CSV_DIALECT)
			for row in reader:
				if row['units'] != '':
					paramDict[row['name']] = json.loads(row['value']) * eval(row['units'])
				else:
					paramDict[row['name']] = json.loads(row['value'])
		setattr(self, attrName, paramDict)
