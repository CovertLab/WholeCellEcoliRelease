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

from wholecell.utils import units

CSV_DIALECT = csv.excel_tab
FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = (
	"compartments.tsv",
	"complexationReactions.tsv",
	"enzymeKinetics.tsv",
	"genes.tsv",
	"metabolites.tsv",
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
	"nutrients.tsv",
	"secretions.tsv",
	"water.tsv",
	"chromosome.tsv",
	"massAtReplicationInitiation.tsv",
	"equilibriumReactions.tsv",
	"foldChanges.tsv",
	"tfIds.tsv",
	"endoRnases.tsv",
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
	)
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = ("parameters.tsv", "mass_parameters.tsv")
CONSTANTS_FILENAME = "constants.tsv"

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
		attrName = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(self, attrName, [])

		with open(file_name, 'rU') as csvfile:
			reader = JsonReader(csvfile, dialect = CSV_DIALECT)
			setattr(self, attrName, [row for row in reader])

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
