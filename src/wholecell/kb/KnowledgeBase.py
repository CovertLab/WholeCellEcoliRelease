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
import Bio.Seq
import Bio.Alphabet.IUPAC


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
		wb = xl.load_workbook(filename = self.dataFileName, use_iterators = True)
		ws = wb.get_sheet_by_name("Genes").iter_rows()

		# Skip the first row
		ws.next()

		self.genes = []
		self.rnas = []
		self.proteins = []
		for row in ws:
			# Gene
			g = {
				"id": row[0].internal_value,
				"name": row[1].internal_value,
				"symbol": row[2].internal_value,
				"type": row[3].internal_value,
				"start": int(row[4].internal_value),
				"len": int(row[5].internal_value),
				"dir": row[6].internal_value == "forward",
				"seq": "",
				"rnaId": row[0].internal_value
			}
			if g["name"] == None:
				g["name"] = ""
			if g["symbol"] == None:
				g["symbol"] = ""
			# import ipdb
			# ipdb.set_trace()
			g["seq"] = self.genomeSeq[(g["start"] - 1) : (g["start"] + g["len"] - 1)]
			if not g["dir"]:
				g["seq"] = Bio.Seq.Seq(g["seq"]).reverse_complement().tostring()
			self.genes.append(g)

			r = {
				"id": g["id"],
				"name": g["name"],
				"type": g["type"],
				"exp": row[7].internal_value,
				"halfLife": row[8].internal_value,
				"seq": "",
				"ntCount": [],
				"mw": -1,
				"geneId": g["id"],
				"monomerId": g["id"] + "_MONOMER"
			}
			# TODO: Figure out why we're taking the complement
			r["seq"] = Bio.Seq.Seq(g["seq"], Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).complement().transcribe().tostring()
			r["ntCount"] = [r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")]
			r["mw"] = 345.20 * r["ntCount"][0] + 321.18 * r["ntCount"][1] + 361.20 * r["ntCount"][2] + 322.17 * r["ntCount"][3] - (len(r["seq"]) - 1) * 17.01
			self.rnas.append(r)

	def loadComplexes(self):
		pass

	def loadReactions(self):
		pass