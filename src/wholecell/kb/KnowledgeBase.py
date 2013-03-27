#!/usr/bin/env python

"""
KnowledgeBase

Whole-cell knowledge base

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/26/2013
"""

import os.path
import openpyxl as xl
import Bio.SeqIO


class KnowledgeBase(object):
	""" KnowledgeBase """

	def __init__(self, dataFileName = "data/KnowledgeBase.xlsx", seqFileName = "data/KnowledgeBase.fna"):
		if not os.path.isfile(dataFileName):
			raise Exception, "%s is missing" % dataFileName
		if not os.path.isfile(seqFileName):
			raise Exception, "%s is missing" % seqFileName

		self.dataFileName = dataFileName
		self.seqFileName = seqFileName

		# Parse data
		self.loadMetabolites()
		self.loadGenome()
		self.loadGenes()
		self.loadComplexes()
		self.loadReactions()

	def loadMetabolites(self):
		wb = xl.load_workbook(filename = self.dataFileName, use_iterators = True)
		ws = wb.get_sheet_by_name("Metabolites").iter_rows()

		# Skip the first row
		ws.next()

		self.metabolites = []
		for row in ws:
			m = {
				"id": row[0].internal_value,
				"name": row[1].internal_value,
				"formula": row[2].internal_value,
				"smiles": row[3].internal_value,
				"charge": row[4].internal_value,
				"mw": row[5].internal_value,
				"hydrophobic": row[6].internal_value == "Y",
				"mediaConc": 0,
				"biomassConc": 0,
				"metabolismNewFlux": 0,
				"metabolismRecyclingFlux": 0,
				"maxExchangeRate": 0
			}
			if m["name"] == None:
				m["name"] = ""
			if m["smiles"] == None:
				m["smiles"] = ""
			if row[7].internal_value != None:
				m["mediaConc"] = row[7].internal_value
			if row[8].internal_value != None:
				m["biomassConc"] = row[8].internal_value
			if row[9].internal_value != None:
				m["metabolismNewFlux"] = row[9].internal_value
			if row[10].internal_value != None:
				m["metabolismRecyclingFlux"] = row[10].internal_value
			if row[11].internal_value != None:
				m["maxExchangeRate"] = row[11].internal_value
			self.metabolites.append(m)

	def loadGenome(self):
		self.translationTable = 4 # E. coli is 11
		self.genomeSeq = Bio.SeqIO.parse(self.seqFileName, "fasta").next().seq.tostring()

	def loadGenes(self):
		pass

	def loadComplexes(self):
		pass

	def loadReactions(self):
		pass