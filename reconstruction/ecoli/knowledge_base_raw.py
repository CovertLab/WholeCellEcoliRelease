"""
KnowledgeBase for Ecoli

Whole-cell knowledge base for Ecoli. Contains all raw, un-fit data processed
directly from CSV flat files.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/11/2015
"""
from __future__ import division

from os import listdir
from os.path import isfile, join
import csv
from reconstruction.spreadsheets import JsonReader

CSV_DIALECT = csv.excel_tab
FLAT_DIR = '/home/users/nruggero/Repos/wcEcoli/reconstruction/ecoli/flat'

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self):
		raw_file_paths = self.find_tsv(FLAT_DIR)
		for file_path in raw_file_paths:
			self.load_tsv(file_path)
		
	def find_tsv(self, file_path):
		return [join(file_path,f) for f in listdir(file_path) if isfile(join(file_path,f)) and f[-4:] == '.tsv']

	def load_tsv(self, file_name):
		attrName = file_name[file_name.rfind('/') + 1 : -4]
		setattr(self, attrName, [])

		with open(file_name) as csvfile:
			print file_name
			reader = JsonReader(csvfile, dialect = CSV_DIALECT)
			for row in reader:
				getattr(self, attrName).append(dict([(x, y) for x,y in row.iteritems()]))