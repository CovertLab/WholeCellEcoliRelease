#!/usr/bin/env python

"""
KnowledgeBase for Ecoli

Whole-cell knowledge base ecoli

@author: Sajia Akhter
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 04/04/2014
"""

import os.path
import csv
import json
import copy
import Bio.SeqIO
import Bio.Seq
import Bio.Alphabet.IUPAC
import Bio.SeqUtils.ProtParam
import re
import numpy


import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, Molecule, Location, Comment, ProteinMonomers, Rna, Metabolite, ProteinComplex, ProteinComplexModified, ProteinMonomerModified, RnaModified, RelationStoichiometry, ProteinComplexReactionRelation,ProteinComplexModifiedReaction,ProteinComplexModReactionRelation, ProteinComplexModReactionEnzyme, ProteinMonomerModifiedReaction, ProteinMonomerModReactionEnzyme, ProteinMonomerModReactionRelation, RnaModifiedReaction, RnaModReactionEnzyme, RnaModifiedReactionRelation, MetaboliteReaction, MetaboliteReactionEnzyme, MetaboliteReactionRelation, MetaboliteBiomass, MetaboliteEquivalentEnzyme, Chromosome, GeneSplices, GeneAbsolutentPosition, EntryPositiveFloatData, GeneType


class KnowledgeBaseEcoli(object):
	""" KnowledgeBase """

	def __init__(self):

		self.dataFileDir = ''
		self.seqFileName = ''

		self.aaWeights = {
			"A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19, "G": 75.07, "H": 155.16, "I": 131.18, "K": 146.19, "L": 131.18,
			"M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09, "T": 119.12, "U": 168.05, "V": 117.15, "W": 204.23,
			"Y": 181.19
    	}
		#TODO: Check casting data type
		# Parse data
		self.loadProducts() # ADDED: for accessing info from other table 
		self.loadComments() # ADDED: for accessing info from other table 
		self.loadCompartments()
		self.loadMetabolites()
		self.loadGenome()
		self.loadGenes()
		self.loadRnas()
		self.loadProteinMonomers() #not dome
		self.createModifiedForms()
		self.loadRelationStoichiometry() # ADDED: for accessing info from other table 
		self.loadComplexes() 
		self.loadReactions()

	def loadProducts(self):

		self.allProducts = {}
		self.allProductType = {} #ADDED for thisType in loadRelationStoichiometry 
	
		#molecule
		all_molecules = Molecule.objects.all()
		if len(all_molecules) <=0:
			raise Exception, "Database Access Error: Cannot access public_molecule table"

		for i in all_molecules:
			self.allProducts[i.id] = i.product 
			self.allProductType[i.product] = '' #updated in RNA, monomer, complex, Metabolite, modifiedForm

	def loadComments(self):

		self.allComments = {}
		
		#molecule
		all_comments = Comment.objects.all()
		if len(all_comments) <=0:
			raise Exception, "Database Access Error: Cannot access public_comment table"

		for i in all_comments:
			self.allComments[i.id] = i.comment_str

		
	def loadCompartments(self):

		self.compartments = []
		self.compIdToAbbrev = {}
		self.dbLocationId = {} # ADDED: for accessing info from other table 

		#location 
		all_locations = Location.objects.all()
		if len(all_locations) <=0:
			raise Exception, "Database Access Error: Cannot access public_location table"

		for i in all_locations:			
			c = {"id": i.location_id, "abbrev": i.abbreviation}

			self.compartments.append(c)
			self.compIdToAbbrev[c["id"]] = c["abbrev"]

			self.dbLocationId[i.id] = c["abbrev"]			


	def loadMetabolites(self):
	
		self.metabolites = []
		
		#biomass		
		biomass = {}
		all_biomass = MetaboliteBiomass.objects.all()	
		if len(all_biomass) <=0:
			raise Exception, "Database Access Error: Cannot access public_MetaboliteBiomass table"
		
		for i in all_biomass:
			temp = {
					"mmol/gDCW":float(i.biomass_concentration),
					"location": self.dbLocationId[i.biomass_location_fk_id]
			}
			if i.metabolite_id_fk_id not in biomass:
				biomass[i.metabolite_id_fk_id] = {'core' : [], 'wildtype' : []}
		
			if i.is_core:
				biomass[i.metabolite_id_fk_id]['core'].append(temp)
			if i.is_wildtype:
				biomass[i.metabolite_id_fk_id]['wildtype'].append(temp)

		#equivalent enz		
		equ_enz = {}
		all_equ_enz = MetaboliteEquivalentEnzyme.objects.all()	
		if len(all_equ_enz) <=0:
			raise Exception, "Database Access Error: Cannot access public_MetaboliteEquivalentEnzyme table"

		for i in all_equ_enz:
			temp = {
					'id' : self.allProducts[i.equivalent_enzyme_id_fk_id], 
					'location' : self.dbLocationId[i.location_fk_id]
			}
			if i.metabolite_id_fk_id not in equ_enz:
				equ_enz[i.metabolite_id_fk_id] = []

			equ_enz[i.metabolite_id_fk_id].append(temp)			

			
		#metabolite
		all_metabolites = Metabolite.objects.all()		
		if len(all_metabolites) <=0:
			raise Exception, "Database Access Error: Cannot access public_Metabolite table"

		for i in all_metabolites:
			
			m = {
					"id": self.allProducts[i.metabolite_id_id].upper(),
					"name": i.name,
					"feistFormula": i.feist_formula,
					"formula7.2": i.ph_formula,
					"charge7.2": int(i.ph_charge),
					"mw7.2": float(i.ph_weight),
					"mediaConc": float(i.media_concentration),
					"biomassInfo": {'core' : [], 'wildtype' : []},
					"maxExchange": float(i.maximum_exchange_rate),
					"fakeMet": i.fake_metabolite,
					"equivEnzIds": [],
					"comments": self.allComments[i.comment_fk_id]
				}
			
			if i.has_biomass:
				m["biomassInfo"] = biomass[i.id]

			if i.has_equivalent_enzyme:
				m["equivEnzIds"] = equ_enz[i.id]		
			
			self.metabolites.append(m)
			self.allProductType[self.allProducts[i.metabolite_id_id]] = 'metabolite' #added

	def loadGenome(self):
		self.translationTable = 11 # E. coli is 11
		
		all_seq = Chromosome.objects.all()
		genome = ''
		for i in all_seq:
			genome = i.sequence
			break
		self.genomeSeq = genome
		

	def loadGenes(self):

		self.genes = []
		self.rnas = []
		self.proteins = []
		self.geneDbIds = {} # ADDED: for rnas and monomers

		#genetype
		genetypes = {}
		all_genetypes = GeneType.objects.all()
		if len(all_genetypes) <=0:
			raise Exception, "Database Access Error: Cannot access public_GeneType table"
		for i in all_genetypes:
			genetypes[i.id] = i.type_gene
		
		#genesplices
		genesplices = {}
		all_genesplices = GeneSplices.objects.all()
		if len(all_genesplices) <=0:
			raise Exception, "Database Access Error: Cannot access public_GeneSplices table"
		for i in all_genesplices:
			genesplices[i.gene_id] = {'start1':int(i.start1),'stop1':int(i.stop1), 											'start2':int(i.start2),'stop2':int(i.stop2)}		

		#EntryPositiveFloatData
		posData = {}
		all_posData = EntryPositiveFloatData.objects.all()
		if len(all_posData) <=0:
			raise Exception, "Database Access Error: Cannot access public_EntryPositiveFloatData table"
		for i in all_posData:
			posData[i.id] = float(i.value)

		#GeneAbsolutentPosition
		genePos = {}
		all_genePos = GeneAbsolutentPosition.objects.all()
		if len(all_genePos) <=0:
			raise Exception, "Database Access Error: Cannot access public_GeneAbsolutentPosition table"
		for i in all_genePos:
			genePos[i.gene_id] = {'pos':int(i.abs_nt_pos),'old':i.old,'new':i.new}

		#Gene
		all_genes = Gene.objects.all()
		if len(all_genes) <=0:
			raise Exception, "Database Access Error: Cannot access public_Gene table"

		for i in all_genes:
			self.geneDbIds[i.id] = i.frame_id # Added for rnas and monomers

			g = {
				"id": i.frame_id,
				"name": str(i.name),
				"symbol": i.symbol,
				"type": genetypes[i.typegene_id],
				"coordinate": int(i.coordinate) - 1, # The coordinates we're given are 1 indexed.
				"length": int(i.length),
				"direction": i.direction,
				"seq": "",
				"rnaId": ""
			}
			g["name"] = g["name"].replace("\\","")
			if g["direction"] == "f":
				g["direction"] = '+'
				g["seq"] = self.genomeSeq[(g["coordinate"]): (g["coordinate"] + g["length"])]
			else:
				g["direction"] = '-'
				g["seq"] = Bio.Seq.Seq(self.genomeSeq[(g["coordinate"] - g["length"] + 1): (g["coordinate"] + 1)]).reverse_complement().tostring()

			if g["type"] == "mRNA":
				g["rnaId"] = g["id"] + "_RNA"
			else:
				g["rnaId"] = self.allProducts[i.productname_id]

			self.genes.append(g)

			# RNA
			r = {
				"id": g["rnaId"],
				"name": "",
				"expression": posData[i.expression_id],
				"modifiedForms": [],
				"unmodifiedForm": None,
				"composition": [],
				"halfLife": posData[i.half_life_id],
				"location": None,
				"seq": "",
				"ntCount": [],
				"mw": -1.0,
				"geneId": g["id"],
				"monomerId": None
			}
			if g["type"] == "mRNA":
				r["name"] = g["name"] + " [RNA]" # else: need to check name in the RNAs file
				r["monomerId"] = self.allProducts[i.productname_id]
				r["location"] = self.compIdToAbbrev["CCO-CYTOSOL"]

				# TODO from DEREK: Uncomment when Nick has fixed json formatting
				# if type(r["halfLife"]) == dict:
				# 	if r["halfLife"]["units"] != "day":
				# 		raise Exception, "Unknown unit!"
				# 	r["halfLife"] = r["halfLife"]["value"] * 24.0 * 60.0 * 60.0

			r["seq"] = Bio.Seq.Seq(g["seq"], Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()
			r["ntCount"] = numpy.array([r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")])
			r["mw"] = 345.20 * r["ntCount"][0] + 321.18 * r["ntCount"][1] + 361.20 * r["ntCount"][2] + 322.17 * r["ntCount"][3] - (len(r["seq"]) - 1) * 17.01
			
			self.rnas.append(r)
			self.allProductType[r["id"]] = 'rna' #added

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
				
				if int(i.splices):
					baseSequence = Bio.Seq.Seq("", Bio.Alphabet.IUPAC.IUPACUnambiguousDNA())
					baseSequence = baseSequence + self.genomeSeq[genesplices[i.id]['start1']-1:genesplices[i.id]['stop1']] + self.genomeSeq[genesplices[i.id]['start2']-1:genesplices[i.id]['stop2']]

					if g["direction"] == "-":
						baseSequence = baseSequence.reverse_complement()
					baseSequence = baseSequence.tostring()

				else:
					baseSequence = g["seq"]

				p["seq"] = Bio.Seq.Seq(baseSequence, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).translate(table = self.translationTable).tostring()

				if int(i.absolute_nt_position):
					pos = genePos[i.id]['pos']
					before = genePos[i.id]['old']
					after = genePos[i.id]['new']
					seqList = list(p["seq"])

					if seqList[pos - 1] != before:
						raise Exception, "Amino acid substitution appears to be incorrect."
					else:
						seqList[pos - 1] = after
					p["seq"] = "".join(seqList)

				p["seq"] = p["seq"][:p["seq"].find('*')]
				tmp = dict([(x, 0) for x in self.aaWeights])

				for aa in tmp: 
					tmp[aa] = p["seq"].count(aa)
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
				self.allProductType[p["id"]] = 'protein' #added		

	def loadRnas(self):

		#RnaModified
		rnamodified = {}
		all_rnamodified = RnaModified.objects.all()
		if len(all_rnamodified) <=0:
			raise Exception, "Database Access Error: Cannot access public_RnaModified table"
		for i in all_rnamodified:
			if i.unmodified_rna_fk_id not in rnamodified:
				rnamodified[i.unmodified_rna_fk_id] = []
			rnamodified[i.unmodified_rna_fk_id].append(str(self.allProducts[i.rna_mod_id]))	

		#rna
		all_rna = Rna.objects.all()
		if len(all_rna) <=0:
			raise Exception, "Database Access Error: Cannot access public_Rna table"

		# rnaId -> location index in self.rnas
		rnaLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self.rnas)])
		
		for i in all_rna:			
			# RNA
			r = {
				"id": self.allProducts[i.frame_id_id],
				"name": i.name,
				"geneId": self.geneDbIds[i.gene_fk_id],
				"location": self.dbLocationId[i.location_fk_id],
				"modifiedForms": [],
				"comments": self.allComments[i.comment_fk_id]
				}
		
			if int(i.is_modified):	
				r["modifiedForms"] = rnamodified[i.id]

			self.rnas[rnaLookup[r["id"]]]["name"] = r["name"]
			self.rnas[rnaLookup[r["id"]]]["location"] = r["location"]
			self.rnas[rnaLookup[r["id"]]]["modifiedForms"] = r["modifiedForms"]
			self.rnas[rnaLookup[r["id"]]]["comments"] = r["comments"]
		

	def loadProteinMonomers(self):

		#ProteinMonomerModified
		monomermod = {}
		all_monomermod = ProteinMonomerModified.objects.all()
		if len(all_monomermod) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinMonomerModified table"
		for i in all_monomermod:
			if i.unmodified_protein_monomer_fk_id not in monomermod:
				monomermod[i.unmodified_protein_monomer_fk_id] = []
			monomermod[i.unmodified_protein_monomer_fk_id].append(str(self.allProducts[i.protein_monomer_mod_id]))

		#ProteinMonomers
		all_monomers = ProteinMonomers.objects.all()
		if len(all_monomers) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinMonomers table"

		# monomerId -> location index in self.proteins
		protLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self.proteins)])

		for i in all_monomers:

			# Monomer
			p = {
				"id": self.allProducts[i.frame_id_id],
				"name": i.name,
				"geneId": self.geneDbIds[i.gene_fk_id],
				"location": self.dbLocationId[i.location_fk_id],
				"modifiedForms": [],
				"comments": self.allComments[i.comment_fk_id]
			}
			if int(i.is_modified):	#TODO: Check after update monomer_modified by Nick
				if i.id in monomermod: 
					p["modifiedForms"] = monomermod[i.id]
				else:
					raise Exception, "modified Monomer Absent %s" % p["id"]
					#print p["id"], i.id, i.name

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
					self.allProductType[rNew["id"]] = 'rna' #added		
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
					self.allProductType[pNew["id"]] = 'protein' #added		
					proteinsToAppend.append(pNew)

		self.proteins.extend(proteinsToAppend)


	def loadRelationStoichiometry(self):

		self.allRelationStoichiometry = {}
		
		all_RelationStoichiometry = RelationStoichiometry.objects.all()
		if len(all_RelationStoichiometry) <=0:
			raise Exception, "Database Access Error: Cannot access public_RelationStoichiometry table"

		for i in all_RelationStoichiometry:
			thisType = self.allProductType[self.allProducts[i.reactant_fk_id]]
			self.allRelationStoichiometry[i.id] = { "coeff": float(i.coefficient), "location": self.dbLocationId[i.location_fk_id], "molecule": self.allProducts[i.reactant_fk_id], "form": "mature", "type":  thisType}


	def loadComplexes(self):
		
		##reaction		
		relation = {}
		all_relation = ProteinComplexReactionRelation.objects.all()
		if len(all_relation) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinComplexReactionRelation table"
		for i in all_relation:
			if i.protein_complex_fk_id not in relation:
				relation[i.protein_complex_fk_id] = []
			relation[i.protein_complex_fk_id].append(i.reactant_relation_id) 

		#proteinComplexesModified
		complexMod = {}
		all_complexMod = ProteinComplexModified.objects.all()
		if len(all_complexMod) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinComplexModified table"
		for i in all_complexMod:
			if i.unmodified_protein_complex_fk_id not in complexMod:
				complexMod[i.unmodified_protein_complex_fk_id] = []
			complexMod[i.unmodified_protein_complex_fk_id].append(str(self.allProducts[i.protein_complex_mod_id]))	

		#proteinComplexes
		all_complex = ProteinComplex.objects.all()
		if len(all_complex) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinComplex table"

		protNew = []		
		for i in all_complex:
			p = {
				"id": self.allProducts[i.protein_complex_id],
				"name": i.name,
				"modifiedForms": [],
				"unmodifiedForm": None,
				"location": self.dbLocationId[i.location_fk_id],
				"composition": [],
				"dir": int(i.reaction_direction),
				"formationProcess": "Complexation",
				"seq": "",
				"aaCount": numpy.zeros(21),
				"ntCount": numpy.zeros(4),
				"mw": -1,
				"geneId": "",
				"rnaId": "",
				"comments": self.allComments[i.comment_fk_id]
			}
			if i.modified_form:	
				p["modifiedForms"] = complexMod[i.id]

			if i.id not in relation:
				raise Exception, "%s protein complex has no reaction" % i.frame_id
			for temp in relation[i.id]:
				t = self.allRelationStoichiometry[temp]
				p["composition"].append(t)


			protNew.append(p)
			self.proteins.append(p)

		self.createModifiedForms()
		
		metDict = dict([(x["id"], x) for x in self.metabolites])
		rnaDict = dict([(x["id"], x) for x in self.rnas])
		protDict = dict([(x["id"], x) for x in self.proteins])
		
		for p in protNew:
			p = [x for x in self.proteins if x["id"] == p["id"]][0]

			for stoichComponent in p["composition"]:
				if stoichComponent["type"] == '': 
					stoichComponent["type"] = self.check_molecule(stoichComponent["molecule"]) #ADDED
				if stoichComponent["molecule"] != p["id"]:
					if stoichComponent["molecule"].upper() in metDict:
						stoichComponent["molecule"] = stoichComponent["molecule"].upper()
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
		
		#self.proteins.extend(protNew)

	def loadReactions(self):

		self.reactions = []
		self.reactionsExchange = []
		
		##		
		relation = {}
		all_relation = MetaboliteReactionRelation.objects.all()
		if len(all_relation) <=0:
			raise Exception, "Database Access Error: Cannot access public_MetaboliteReactionRelation table"
		for i in all_relation:
			if i.metabolite_reaction_fk_id not in relation:
				relation[i.metabolite_reaction_fk_id] = []
			relation[i.metabolite_reaction_fk_id].append(i.reactant_relation_id) 
		
		##
		enz = {}
		all_enz = MetaboliteReactionEnzyme.objects.all()
		if len(all_enz) <=0:
			raise Exception, "Database Access Error: Cannot access public_MetaboliteReactionEnzyme table"
		for i in all_enz:
			if i.metabolite_reaction_fk_id not in enz:
				enz[i.metabolite_reaction_fk_id] = []
			enz[i.metabolite_reaction_fk_id].append(str(self.allProducts[i.enzyme_fk_id])) 
		
		##
		all_metaboliteReaction = MetaboliteReaction.objects.all()
		if len(all_metaboliteReaction) <=0:
			raise Exception, "Database Access Error: Cannot access public_metaboliteReaction table"
		

		for i in all_metaboliteReaction:
			r = {
					"id": i.frame_id,
					"name": i.name,
					"process": "Metabolism",
					"ec": i.ec,
					"dir": int(i.reaction_direction),
					"stoichiometry": [],
					"catBy": [],
					"ub": float(i.upper_bound),
					"lb": float(i.lower_bound)
				}

			if i.id in enz:
				r["catBy"] = enz[i.id]
			else:
				r["catBy"] = None

			if r["name"] == None: r["name"] = ""
			if r["ec"] == None: r["ec"] = ""
			
			if i.id not in relation:
				raise Exception, "%s Metabolite has no reaction" % i.frame_id
			for temp in relation[i.id]:
				t = self.allRelationStoichiometry[temp]
				if t["type"] == '':
					t["type"] = self.check_molecule(t["molecule"])
				t["molecule"] = t["molecule"].upper()
				r["stoichiometry"].append(t)
	
			protList = [p["id"] for p in self.proteins]
			if r["catBy"] != None:
				for temp in r["catBy"]:
					if temp not in protList:
						raise Exception, "Undefined protein: %s." % temp
			
			self.reactions.append(r)


	def check_molecule(self, mol):
		thisType = ""
		'''		
		if any(x["id"] == mol.upper() for x in self.metabolites):
			thisType = "metabolite"
		elif any(x["id"] == mol for x in self.rnas):
			thisType = "rna"
		el
		'''
		if any(x["id"] == mol for x in self.proteins):
			thisType = "protein"
		else:
			raise Exception, "Undefined molecule: %s." % (mol)
		return thisType
				
				
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

