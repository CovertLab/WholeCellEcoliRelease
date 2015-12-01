"""
ValidationDataRawEcoli

Raw validation data for Ecoli. Contains raw data processed
directly from TSV flat files.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/30/2015
"""

import os

FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = (
	"taniguichi2010_table_6.tsv",
	"houser2015_javier_table.tsv",
	# os.path.join("rna_seq_data","rnaseq_seal_rpkm_std.tsv"),
	)

class ValidationDataRawEcoli(object):
	""" ValidationDataEcoli """

	def __init__(self):
		# Load raw data from TSV files
		for filename in LIST_OF_DICT_FILENAMES:
			self._load_tsv(os.path.join(FLAT_DIR, filename))

	def _load_tsv(self, file_name):
		attrName = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(self, attrName, [])


		with open(file_name, 'r') as csvfile:
			rows = []
			for idx, row in enumerate(csvfile):
				if idx == 0:
					header_names = [x.strip("\n") for x in row.split("\t")]
					continue
				columns = [x.strip("\n") for x in row.split("\t")]
				outputrow = dict(zip(header_names,columns))
				rows.append(outputrow)
			setattr(self, attrName, rows)