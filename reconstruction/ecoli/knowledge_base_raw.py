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

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self):
		pass

	def load_tsv(self, file_name):
		# Nick: you probably don't need to make this a method, since the
		# JsonReader class handles most of the parsing.  But the pattern below
		# should work. - John
		with open(file_name) as csvfile:
			reader = JsonReader(csvfile, dialect = csv.excel_tab)
