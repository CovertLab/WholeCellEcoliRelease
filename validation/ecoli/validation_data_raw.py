"""
ValidationDataRawEcoli

Raw validation data for Ecoli. Contains raw data processed
directly from TSV flat files.
"""

from __future__ import absolute_import, division, print_function

import os
from reconstruction.spreadsheets import read_tsv

FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = (
	"taniguichi2010_table_6.tsv",
	"houser2015_javier_table.tsv",
	"wisniewski2014_supp2.tsv",
	"schmidt2015_javier_table.tsv",
	"toya_2010_central_carbon_fluxes.tsv",
	"dna_footprint_sizes.tsv",
	"essential_genes.tsv",
	"geneFunctions.tsv",
	)

class ValidationDataRawEcoli(object):
	""" ValidationDataRawEcoli """

	def __init__(self):
		"""Load raw data from TSV files"""
		for filename in LIST_OF_DICT_FILENAMES:
			self._load_tsv(os.path.join(FLAT_DIR, filename))

	def _load_tsv(self, file_name):
		attrName = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(self, attrName, [])

		rows = read_tsv(file_name)
		setattr(self, attrName, rows)
