"""
ValidationData for Ecoli

Whole-cell knowledge base for Ecoli, valdiation data only. Contains all raw, un-fit data processed
directly from CSV flat files.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/30/2015
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
	"taniguichi2010_table_6.tsv",
	# os.path.join("rna_seq_data","rnaseq_seal_rpkm_std.tsv"),
	)
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = ("parameters.tsv", "mass_parameters.tsv")
CONSTANTS_FILENAME = "constants.tsv"

class ValidationDataEcoli(object):
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
