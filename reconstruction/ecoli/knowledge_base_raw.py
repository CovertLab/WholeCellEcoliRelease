"""
KnowledgeBase for Ecoli
Whole-cell knowledge base for Ecoli. Contains all raw, un-fit data processed
directly from CSV flat files.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""
from __future__ import absolute_import, division, print_function

import os
import csv
from reconstruction.spreadsheets import JsonReader
import json
from itertools import ifilter

from wholecell.utils import units  # used by eval()

CSV_DIALECT = csv.excel_tab
FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = (
	"biomass.tsv",
	"compartments.tsv",
	"complexationReactions.tsv",
	"disabledKineticReactions.tsv",
	"dryMassComposition.tsv",
	"endoRnases.tsv",
	"enzymeKinetics.tsv",
	"equilibriumReactions.tsv",
	"foldChanges.tsv",
	"full_chromosome.tsv",
	"genes.tsv",
	"growthRateDependentParameters.tsv",
	"massAtReplicationInitiation.tsv",
	"metabolites.tsv",
	"metaboliteConcentrations.tsv",
	"modificationReactions.tsv",
	"modifiedForms.tsv",
	"modifiedFormsStoichiometry.tsv",
	"modifiedRnas.tsv",
	"polymerized.tsv",
	"ppgpp_fc.tsv",
	"ppgpp_regulation.tsv",
	"previousBiomassFluxes.tsv",
	"promoters.tsv",
	"protein_half_lives.tsv",
	"proteinComplexes.tsv",
	"proteins.tsv",
	"reactions.tsv",
	"rnas.tsv",
	"secretions.tsv",
	"sequence_motifs.tsv",
	"terminators.tsv",
	"tfIds.tsv",
	"tfOneComponentBound.tsv",
	"transcriptionUnits.tsv",
	"translationEfficiency.tsv",
	"transport_reactions.tsv",
	"twoComponentSystemTemplates.tsv",
	"twoComponentSystems.tsv",
	"water.tsv",
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
	os.path.join("common_names", "mrna.tsv"),
	os.path.join("common_names", "polypeptides.tsv"),
	os.path.join("common_names", "protein-complexes.tsv"),
	os.path.join("common_names", "reactions.tsv"),
	os.path.join("common_names", "rnas.tsv"),
	)
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = ("parameters.tsv", "mass_parameters.tsv")

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
		for subPath in file_name[len(dir_name) + 1 : ].split(os.path.sep)[:-1]:
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
				value = json.loads(row['value'])
				if row['units'] != '':
					# `eval()` the units [risky!] then strip it to just a unit
					# since `a_list * a_float` (like `1.0 [1/s]`) fails, and
					# `a_list * an_int` repeats the list, which is also broken.
					unit = eval(row['units'])   # risky!
					unit = units.getUnit(unit)  # strip
					value = value * unit
				paramDict[row['name']] = value

		setattr(self, attrName, paramDict)
