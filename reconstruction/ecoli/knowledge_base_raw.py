
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
	"biomass.tsv",
	"chromosome.tsv",
	"compartments.tsv",
	"complexationReactions.tsv",
	"disabledKineticReactions.tsv",
	"dryMassComposition.tsv",
	"dryMassComposition_alternateProtein.tsv",
	"dryMassComposition_alternateProteinAndRna.tsv",
	"dryMassComposition_alternateRna.tsv",
	"endoRnases.tsv",
	"enzymeKinetics.tsv",
	"equilibriumReactions.tsv",
	"foldChanges.tsv",
	"genes.tsv",
	"growthRateDependentParameters.tsv",
	"growthRateDependentParameters_alternateRibosomeActivity.tsv",
	"massAtReplicationInitiation.tsv",
	"metabolites.tsv",
	"metaboliteConcentrations.tsv",
	"modificationReactions.tsv",
	"modifiedForms.tsv",
	"modifiedFormsStoichiometry.tsv",
	"modifiedRnas.tsv",
	"polymerized.tsv",
	"previousBiomassFluxes.tsv",
	"promoters.tsv",
	"protein_half_lives.tsv",
	"proteinComplexes.tsv",
	"proteins.tsv",
	"reactions.tsv",
	"ribosomal_protein_transcripts.tsv",
	"rnas.tsv",
	"rnas_alternate_half_lives_without_kas.tsv",
	"secretions.tsv",
	"terminators.tsv",
	"tfIds.tsv",
	"tfOneComponentBound.tsv",
	"transcriptionUnits.tsv",
	"translationEfficiency.tsv",
	"translationEfficiency_alternate.tsv",
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
	os.path.join("rna_seq_data","alternate_rna_seq.tsv"),
	os.path.join("condition", "tf_condition.tsv"),
	os.path.join("condition", "condition_defs.tsv"),
	os.path.join("condition", "nutrient", "minimal.tsv"),
	os.path.join("condition", "nutrient", "minimal_acetate.tsv"),
	os.path.join("condition", "nutrient", "minimal_fumarate.tsv"),
	os.path.join("condition", "nutrient", "minimal_malate.tsv"),
	os.path.join("condition", "nutrient", "minimal_minus_calcium.tsv"),
	os.path.join("condition", "nutrient", "minimal_minus_magnesium.tsv"),
	os.path.join("condition", "nutrient", "minimal_minus_oxygen.tsv"),
	os.path.join("condition", "nutrient", "minimal_minus_phosphate.tsv"),
	os.path.join("condition", "nutrient", "minimal_no_glucose.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_amino_acids.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_cytidine.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_ferric.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_gallate.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_indole.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_nitrate.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_nitrite.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_quercetin.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_tungstate.tsv"),
	os.path.join("condition", "nutrient", "minimal_succinate.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_sam.tsv"),
	os.path.join("condition", "nutrient", "minimal_plus_arabinose.tsv"),
	os.path.join("condition", "timeseries", "000000_basal.tsv"),
	os.path.join("condition", "timeseries", "000001_cut_glucose.tsv"),
	os.path.join("condition", "timeseries", "000002_add_aa.tsv"),
	os.path.join("condition", "timeseries", "000003_aa.tsv"),
	os.path.join("condition", "timeseries", "000004_oxygen_absent.tsv"),
	os.path.join("condition", "timeseries", "000005_indole_present.tsv"),
	os.path.join("condition", "timeseries", "000006_tungstate_present.tsv"),
	os.path.join("condition", "timeseries", "000007_quercetin_present.tsv"),
	os.path.join("condition", "timeseries", "000008_gallate_present.tsv"),
	os.path.join("condition", "timeseries", "000009_succinate_carbon_source.tsv"),
	os.path.join("condition", "timeseries", "000010_acetate_carbon_source.tsv"),
	os.path.join("condition", "timeseries", "000011_fumarate_carbon_source.tsv"),
	os.path.join("condition", "timeseries", "000012_malate_carbon_source.tsv"),
	os.path.join("condition", "timeseries", "000013_nitrate_present.tsv"),
	os.path.join("condition", "timeseries", "000014_nitrite_present.tsv"),
	os.path.join("condition", "timeseries", "000015_calcium_absent.tsv"),
	os.path.join("condition", "timeseries", "000016_magnesium_absent.tsv"),
	os.path.join("condition", "timeseries", "000017_phosphate_absent.tsv"),
	os.path.join("condition", "timeseries", "000018_cut_oxygen.tsv"),
	os.path.join("condition", "timeseries", "000019_add_arabinose.tsv"),
	os.path.join("condition", "timeseries", "000020_add_oxygen.tsv"),
	os.path.join("condition", "timeseries", "000021_add_indole.tsv"),
	os.path.join("condition", "timeseries", "000022_cut_indole.tsv"),
	os.path.join("condition", "timeseries", "000023_add_calcium.tsv"),
	os.path.join("condition", "timeseries", "000024_cut_calcium.tsv"),
	os.path.join("condition", "timeseries", "000025_cut_aa.tsv"),
	os.path.join("condition", "timeseries", "000026_add_and_cut_aa.tsv"),
	)
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = (
	"parameters.tsv", "mass_parameters.tsv", "mass_parameters_alternate.tsv")
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
