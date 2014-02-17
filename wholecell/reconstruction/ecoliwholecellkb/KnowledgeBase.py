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
import csv
import json
import copy
import Bio.SeqIO
import Bio.Seq
import Bio.Alphabet.IUPAC
import Bio.SeqUtils.ProtParam
import re
import numpy

class KnowledgeBase(object):
	""" KnowledgeBase """

	def __init__(self, dataFileDir = "/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli", seqFileName = "/home/wholecell/ecoliwholecellkb/public/fixtures/sequence.fasta"):
		if not os.path.isdir(dataFileDir):
			raise Exception, "%s is missing." % dataFileName
		if not os.path.isfile(seqFileName):
			raise Exception, "%s is missing." % seqFileName

		self.dataFileDir = dataFileDir
		self.seqFileName = seqFileName

		# Borrowed from BioPython
		# They have the Selenocysteine (U) value commented out in IUPACData.py, so we can't use their functions
		# Fortunately they are simple functions with source code available at: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-pysrc.html
		self.aaWeights = {
			"A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19, "G": 75.07, "H": 155.16, "I": 131.18, "K": 146.19, "L": 131.18,
			"M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09, "T": 119.12, "U": 168.05, "V": 117.15, "W": 204.23,
			"Y": 181.19
    	}

		# Parse data
		self.loadCompartments()
		self.loadMetabolites()
		self.loadGenome()
		self.loadGenes()
		self.loadRnas()
		self.loadProteinMonomers()
		self.createModifiedForms()
		self.loadComplexes()
		self.loadReactions()

	def loadCompartments(self):
		fileName = self.dataFileDir + os.sep + "locations.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing" % fileName

		self.compartments = []
		self.compIdToAbbrev = {}
		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "abbrev"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")

			for row in dr:
				c = {
					"id": row["id"],
					"abbrev": row["abbrev"]
				}
				self.compartments.append(c)
				self.compIdToAbbrev[c["id"]] = c["abbrev"]

	def loadMetabolites(self):
		fileName = self.dataFileDir + os.sep + "metabolites.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing" % fileName

		self.metabolites = []
		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "name", "formula7.2", "charge7.2", "mw7.2", "feistFormula", "mediaConc", "biomassInfo", "maxExchange", "fakeMet", "equivEnzIds", "comments"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")

			for row in dr:
				# Metabolite
				m = {
					"id": row["id"].upper(),
					"name": row["name"],
					"feistFormula": row["feistFormula"],
					"formula7.2": row["formula7.2"],
					"charge7.2": int(row["charge7.2"]),
					"mw7.2": float(row["mw7.2"]),
					"mediaConc": 0.,
					"biomassInfo": {'core' : [], 'wildtype' : []},
					"maxExchange": 0.,
					"fakeMet": False,
					"equivEnzIds": [],
					"comments": ""
				}
				if row["mediaConc"]: m["mediaConc"] = float(row["mediaConc"])

				biomassInformation = json.loads(row["biomassInfo"])
				if biomassInformation["core"] != []:
					m["biomassInfo"]["core"] = biomassInformation["core"]
				if biomassInformation["wildtype"] != []:
					m["biomassInfo"]["wildtype"] = biomassInformation["wildtype"]
				if row["maxExchange"]: m["maxExchange"] = float(row["maxExchange"])
				if row["fakeMet"]:
					if row['fakeMet'] == 'True':
						m["fakeMet"] = True
					else:
						m['fakeMet'] = False
				if row["equivEnzIds"] != "[]":
					rawEquivEnzIds = json.loads(row["equivEnzIds"])
					for enz in rawEquivEnzIds:
						m["equivEnzIds"].append({'id' : enz[:-3], 'location' : enz[-2:-1]})
				if row["comments"]: m["comments"] = row["comments"].strip()

				self.metabolites.append(m)

	def loadGenome(self):
		self.translationTable = 11 # E. coli is 11
		self.genomeSeq = Bio.SeqIO.parse(self.seqFileName, "fasta").next().seq.tostring()

	def loadGenes(self):
		fileName = self.dataFileDir + os.sep + "genes.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing" % fileName

		self.genes = []
		self.rnas = []
		self.proteins = []
		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "name", "symbol", "type", "coordinate", "length", "direction", "expression", "halfLife", "product", "spliceInfo", "substInfo", "comments"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")

			for row in dr:
				# Gene
				g = {
					"id": row["id"],
					"name": row["name"],
					"symbol": row["symbol"],
					"type": row["type"],
					"coordinate": int(row["coordinate"]) - 1, # The coordinates we're given are 1 indexed.
					"length": int(row["length"]),
					"direction": row["direction"],
					"seq": "",
					"rnaId": ""
				}

				if g["direction"] == "+":
					g["seq"] = self.genomeSeq[(g["coordinate"]): (g["coordinate"] + g["length"])]
				else:
					g["seq"] = Bio.Seq.Seq(self.genomeSeq[(g["coordinate"] - g["length"] + 1): (g["coordinate"] + 1)]).reverse_complement().tostring()

				if g["type"] == "mRNA":
					g["rnaId"] = g["id"] + "_RNA"
				else:
					g["rnaId"] = row["product"]

				self.genes.append(g)

				# RNA
				r = {
					"id": g["rnaId"],
					"name": "",
					"expression": float(row["expression"]),
					"modifiedForms": [],
					"unmodifiedForm": None,
					"composition": [],
					"halfLife": -1.0,
					"location": None,
					"seq": "",
					"ntCount": [],
					"mw": -1.0,
					"geneId": g["id"],
					"monomerId": None
				}
				if g["type"] == "mRNA":
					r["name"] = g["name"] + " [RNA]" # else: need to check name in the RNAs file
					r["monomerId"] = row["product"]
					r["location"] = self.compIdToAbbrev["CCO-CYTOSOL"]
		
				# TODO: Uncomment when Nick has fixed json formatting
				r["halfLife"] = float(json.loads(row["halfLife"]))
				# if type(r["halfLife"]) == dict:
				# 	if r["halfLife"]["units"] != "day":
				# 		raise Exception, "Unknown unit!"
				# 	r["halfLife"] = r["halfLife"]["value"] * 24.0 * 60.0 * 60.0

				r["seq"] = Bio.Seq.Seq(g["seq"], Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()
				r["ntCount"] = numpy.array([r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")])
				r["mw"] = 345.20 * r["ntCount"][0] + 321.18 * r["ntCount"][1] + 361.20 * r["ntCount"][2] + 322.17 * r["ntCount"][3] - (len(r["seq"]) - 1) * 17.01
				self.rnas.append(r)

				if g["type"] == "mRNA":
					p = {
						"id": r["monomerId"],
						"name": "", # Get from monomers file
						"modifiedForms": [],
						"unmodifiedForm": None,
						"location": None,
						"composition": [],
						"formationProcess": "",
						"seq": "",
						"aaCount": numpy.zeros(21),
						"ntCount": numpy.zeros(4),
						"mw": -1,
						"geneId": g["id"],
						"rnaId": g["rnaId"]
					}
					if row["spliceInfo"] != "[]":
						baseSequence = Bio.Seq.Seq("", Bio.Alphabet.IUPAC.IUPACUnambiguousDNA())
						for splice in json.loads(row["spliceInfo"]):
							baseSequence += self.genomeSeq[(splice[0] - 1): splice[1]]

						if g["direction"] == "-":
							baseSequence = baseSequence.reverse_complement()

						baseSequence = baseSequence.tostring()
					else:
						baseSequence = g["seq"]
					p["seq"] = Bio.Seq.Seq(baseSequence, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).translate(table = self.translationTable).tostring()

					if row["substInfo"] != "[]":
						pos, before, after = json.loads(row["substInfo"])
						seqList = list(p["seq"])

						if seqList[pos - 1] != before:
							raise Exception, "Amino acid substitution appears to be incorrect."
						else:
							seqList[pos - 1] = after
						p["seq"] = "".join(seqList)

					p["seq"] = p["seq"][:p["seq"].find('*')]

					tmp = dict([(x, 0) for x in self.aaWeights])
					for aa in tmp: tmp[aa] = p["seq"].count(aa)
					p["aaCount"] = numpy.array([tmp["A"], tmp["R"], tmp["N"], tmp["D"], tmp["C"],
									tmp["E"], tmp["Q"], tmp["G"], tmp["H"], tmp["I"],
									tmp["L"], tmp["K"], tmp["M"], tmp["F"], tmp["P"],
									tmp["U"], tmp["S"], tmp["T"], tmp["W"], tmp["Y"], tmp["V"]
									])

					water = 18.02
					aaWeights = {}
					for k in self.aaWeights: aaWeights[k] = self.aaWeights[k] - water
					p["mw"] = water
					for aa in p["seq"]: p["mw"] += aaWeights[aa]

					self.proteins.append(p)

	def loadRnas(self):
		fileName = self.dataFileDir + os.sep + "rna.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing." % fileName

		# rnaId -> location index in self.rnas
		rnaLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self.rnas)])

		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "name", "geneId", "location", "modifiedForms", "comments"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")
			
			for row in dr:
				# RNA
				r = {
					"id": row["id"],
					"name": row["name"],
					"geneId": row["geneId"],
					"location": self.compIdToAbbrev[json.loads(row["location"])[0]],
					"modifiedForms": json.loads(row["modifiedForms"]),
					"comments": row["comments"]
				}

				self.rnas[rnaLookup[r["id"]]]["name"] = r["name"]
				self.rnas[rnaLookup[r["id"]]]["location"] = r["location"]
				self.rnas[rnaLookup[r["id"]]]["modifiedForms"] = r["modifiedForms"]
				self.rnas[rnaLookup[r["id"]]]["comments"] = r["comments"]

	def loadProteinMonomers(self):
		fileName = self.dataFileDir + os.sep + "proteinMonomers.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing." % fileName

		# monomerId -> location index in self.proteins
		protLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self.proteins)])

		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "name", "geneId", "location", "modifiedForms", "comments"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")
			
			for row in dr:
				# Monomer
				p = {
					"id": row["id"],
					"name": row["name"],
					"geneId": row["geneId"],
					"location": self.compIdToAbbrev[json.loads(row["location"])[0]],
					"modifiedForms": json.loads(row["modifiedForms"]),
					"comments": row["comments"]
				}

				self.proteins[protLookup[p["id"]]]["name"] = p["name"]
				self.proteins[protLookup[p["id"]]]["location"] = p["location"]
				self.proteins[protLookup[p["id"]]]["modifiedForms"] = p["modifiedForms"]
				self.proteins[protLookup[p["id"]]]["comments"] = p["comments"]

	def createModifiedForms(self):
		rnaIds = [x["id"] for x in self.rnas]
		rnasToAppend = []
		for r in self.rnas:
			for modForm in r["modifiedForms"]:
				if modForm not in rnaIds:	# Do this check so that we can call the function multiple times and not re-create entries
					rNew = dict(r)
					rNew["id"] = modForm
					rNew["modifiedForms"] = []
					rNew["unmodifiedForm"] = r["id"]
					rNew["composition"] = []
					rNew["mw"] = -1.0 	# TODO: Need to get this
					rNew["expression"] = 0.
					rnasToAppend.append(rNew)

		self.rnas.extend(rnasToAppend)

		protIds = [x["id"] for x in self.proteins]
		proteinsToAppend = []
		for p in self.proteins:
			for modForm in p["modifiedForms"]:
				if modForm not in protIds:	# Do this check so that we can call the function multiple times and not re-create entries
					pNew = dict(p)
					pNew["id"] = modForm
					pNew["modifiedForms"] = []
					pNew["unmodifiedForm"] = p["id"]
					pNew["composition"] = []
					pNew["mw"] = -1.0 	# TODO: Need to get this
					proteinsToAppend.append(pNew)

		self.proteins.extend(proteinsToAppend)

	def loadComplexes(self):
		fileName = self.dataFileDir + os.sep + "proteinComplexes.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing." % fileName

		protNew = []
		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "name", "location", "compositionStr", "compositionDict", "modifiedForms", "formationProcess", "comments"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")
			
			for row in dr:
				p = {
					"id": row["id"],
					"name": row["name"],
					"modifiedForms": json.loads(row["modifiedForms"]),
					"unmodifiedForm": None,
					"location": self.compIdToAbbrev[json.loads(row["location"])[0]],
					"composition": row["compositionStr"],
					"formationProcess": row["formationProcess"],
					"seq": "",
					"aaCount": numpy.zeros(21),
					"ntCount": numpy.zeros(4),
					"mw": -1,
					"geneId": "",
					"rnaId": "",
					"comments": ""
				}
				protNew.append(p)
				self.proteins.append(p)

		self.createModifiedForms()

		metDict = dict([(x["id"], x) for x in self.metabolites])
		rnaDict = dict([(x["id"], x) for x in self.rnas])
		protDict = dict([(x["id"], x) for x in self.proteins])

		for p in protNew:
			p = [x for x in self.proteins if x["id"] == p["id"]][0]

			p["composition"] = self.parseReaction(p["composition"])[0]
			for stoichComponent in p["composition"]:
				if stoichComponent["molecule"] != p["id"]:
					if stoichComponent["molecule"] in metDict:
						subunitMw = metDict[stoichComponent["molecule"]]["mw7.2"]
					elif stoichComponent["molecule"] in rnaDict:
						subunitMw = rnaDict[stoichComponent["molecule"]]["mw"]
						p["ntCount"] -= stoichComponent["coeff"] * rnaDict[stoichComponent["molecule"]]["ntCount"]
					elif stoichComponent["molecule"] in protDict:
						subunitMw = protDict[stoichComponent["molecule"]]["mw"]
						p["aaCount"] -= stoichComponent["coeff"] * protDict[stoichComponent["molecule"]]["aaCount"]
					else:
						raise Exception, "Undefined subunit: %s." % stoichComponent["molecule"]

					p["mw"] -= stoichComponent["coeff"] * subunitMw

#		self.proteins.extend(protNew)

	def loadReactions(self):
		fileName = self.dataFileDir + os.sep + "reactions.csv"
		if not os.path.isfile(fileName):
			raise Exception, "%s is missing." % fileName

		self.reactions = []
		self.reactionsExchange = []
		with open(fileName, "r") as csvfile:

			# Skip the first row
			csvfile.next()

			fieldnames = ["id", "name", "process", "ec", "stoichiometry", "enzyme", "ub", "lb", "comments"]
			dr = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = "\t")
			
			for row in dr:
				if row["id"][:9] == "FEIST_EX_" or row["id"][:9] == "FEIST_DM_":
					parseReaction = self.parseLeakReaction
				else:
					parseReaction = self.parseReaction
				r = {
					"id": row["id"],
					"name": row["name"],
					"process": row["process"],
					"ec": row["ec"],
					"dir": "",
					"stoichiometry": [],
					"catBy": json.loads(row["enzyme"]),
					"ub": float(row["ub"]),
					"lb": float(row["lb"])
				}

				if r["name"] == None: r["name"] = ""
				if r["ec"] == None: r["ec"] = ""
				r["stoichiometry"], r["dir"] = parseReaction(row["stoichiometry"])
				if r["catBy"] != None: r["catBy"] = [str(x) for x in r["catBy"]]
				
				protList = [p["id"] for p in self.proteins]
				if r["catBy"] != None:
					for enz in r["catBy"]:
						if enz not in protList:
							raise Exception, "Undefined protein: %s." % row["enzyme"]

				self.reactions.append(r)

	def parseReaction(self, reactionStr):
		match = re.match("^\[(?P<comp>.*?)\][ ]{0,1}: (?P<stoich>.*)$", reactionStr)
		if match != None:
			globalComp = match.group("comp")
			stoich = match.group("stoich")
		else:
			globalComp = ""
			stoich = reactionStr

		match = re.match("^(?P<lefts>.*) (?P<dir><*((==)|(--))>*) (?P<rights>.*)$", stoich)
		if match == None:
			raise Exception, "Invalid stoichiometry: %s." % (stoich)

		if match.group("dir") == "==>" or match.group("dir") == "-->":
			reactionDir = 1
		elif match.group("dir") == "<==" or match.group("dir") == "<--":
			reactionDir = -1
		elif match.group("dir") == "<==>" or match.group("dir") == "<-->":
			reactionDir = 0

		stoich = []

		lefts = match.group("lefts").split(" + ")
		for componentStr in lefts:
			coeff, mol, form, comp, thisType = self.parseReactionComponent(componentStr, globalComp)
			stoich.append({ "coeff": -coeff, "location": comp, "molecule": mol, "form": form, "type": thisType })

		rights = match.group("rights").split(" + ")
		for componentStr in rights:
			coeff, mol, form, comp, thisType = self.parseReactionComponent(componentStr, globalComp)
			stoich.append({ "coeff": coeff, "location": comp, "molecule": mol, "form": form, "type": thisType })

		return stoich, reactionDir

	def parseLeakReaction(self, reactionStr):
		match = re.match("^\[(?P<comp>.*?)\][ ]{0,1}: (?P<stoich>.*)$", reactionStr)
		if match != None:
			globalComp = match.group("comp")
			stoich = match.group("stoich")
		else:
			globalComp = ""
			stoich = reactionStr

		match = re.match("^(?P<lefts>.*?)[ ]{0,1}((?P<dir><*((==)|(--))>*))*$", stoich)
		if match == None:
			raise Exception, "Invalid stoichiometry: %s." % (stoich)

		if match.group("dir") == "==>" or match.group("dir") == "-->":
			reactionDir = 1
		elif match.group("dir") == "<==" or match.group("dir") == "<--":
			raise Exception, "Invalid stoichiometry for a leak reaction: %s." % (stoich)
		elif match.group("dir") == "<==>" or match.group("dir") == "<-->":
			reactionDir = 0
		elif match.group("dir") == None:
			reactionDir = 0

		stoich = []

		lefts = match.group("lefts").split(" + ")
		for componentStr in lefts:
			coeff, mol, form, comp, thisType = self.parseReactionComponent(componentStr, globalComp)
			stoich.append({ "coeff": -coeff, "location": comp, "molecule": mol, "form": form, "type": thisType })

		if len(stoich) > 1:
			raise Exception, "Leak reaction can only have one metabolite: %s." % (stoich)

		return stoich, reactionDir

	def parseReactionComponent(self, componentStr, globalComp):
		if globalComp == "":
			tmp = re.match("^(?P<coeff>\(\d*\.*\d*\) )*(?P<mol>.+?)(?P<form>:.+?)*(\[(?P<comp>.+)\])*$", componentStr)
			if tmp == None:
				raise Exception, "Invalid stoichiometry: %s." % (componentStr)
			if tmp.group("coeff") == None:
				coeff = 1.0
			else:
				coeff = float(tmp.group("coeff")[1:-2])

			mol = tmp.group("mol")

			if tmp.group("form") == None:
				form = "mature"
			else:
				form = tmp.group("form")[1:]

			if tmp.group("comp") == None:
				comp = "c"
			else:
				comp = tmp.group("comp")
		else:
			tmp = re.match("^(?P<coeff>\(\d*\.*\d*\) )*(?P<mol>.+?)(?P<form>:.+?)*(\[(?P<comp>.+)\])*$", componentStr)
			if tmp == None:
				raise Exception, "Invalid stoichiometry: %s." % (componentStr)
			if tmp.group("coeff") == None:
				coeff = 1.0
			else:
				coeff = float(tmp.group("coeff")[1:-2])

			mol = tmp.group("mol")

			if tmp.group("form") == None:
				form = "mature"
			else:
				form = tmp.group("form")[1:]

			if tmp.group("comp") != None:
				if tmp.group("comp") != globalComp:
					raise Exception, "Metabolite compartment %s does not match global compartment %s." % (tmp.group("comp"), globalComp)

			comp = globalComp
		
		if any(x["id"] == mol.upper() for x in self.metabolites):
			mol = mol.upper()
			thisType = "metabolite"
		elif any(x["id"] == mol for x in self.rnas):
			thisType = "rna"
		elif any(x["id"] == mol for x in self.proteins):
			thisType = "protein"
		else:
			raise Exception, "Undefined molecule: %s." % (mol)
		
		return coeff, mol, form, comp, thisType

	def calcKCat(self, enzId, vMax, units):
		if enzId == None or vMax == None:
			return numpy.NaN

		if units == "U/mg":
			prot = next((x for x in self.proteins if x["id"] == enzId), None)
			if prot == None:
				raise Exception, "Undefined enzyme: %s." % (enzId)
			return vMax / 60.0 * 1e-3 * prot["mw"]
		elif units == "1/min":
			return vMax / 60.0
		else:
			raise Exception, "Invalid kCat units: %s." % (units)
