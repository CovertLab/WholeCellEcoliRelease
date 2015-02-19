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

CSV_DIALECT = csv.excel_tab
FLAT_DIR = os.path.join("reconstruction", "ecoli", "flat")
FILENAMES = (
	"compartments.tsv",
	"complexationReactions.tsv",
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
	)
SEQUENCE_FILE = 'sequence.fasta'

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self):
		# Load raw data from TSV files
		for filename in FILENAMES:
			self._load_tsv(os.path.join(FLAT_DIR, filename))

		self.genome_sequence = self._load_sequence(os.path.join(FLAT_DIR, SEQUENCE_FILE))


	def _load_tsv(self, file_name):
		attrName = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(self, attrName, [])

		with open(file_name) as csvfile:
			reader = JsonReader(csvfile, dialect = CSV_DIALECT)
			for row in reader:
				getattr(self, attrName).append(dict([(x, y) for x,y in row.iteritems()]))

	def _load_sequence(self, file_path):
		from Bio import SeqIO
		with open(file_path, "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				return record.seq
