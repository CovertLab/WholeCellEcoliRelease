#!/usr/bin/env python

"""
KnowledgeBase for Ecoli

Whole-cell knowledge base ecoli

@author: Sajia Akhter
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 04/04/2014
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/14/2014
"""
from __future__ import division
import numpy # TODO: change to import numpy as np
import collections

import os
import sys

# Set Django environmental variable
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb_project.ecoliwholecellkb.settings'

import wholecell.utils.config
sys.path.append(str(os.path.expanduser(wholecell.utils.config.KNOWLEDGEBASE_PACKAGE_DIR)))
import ecoliwholecellkb_project.ecoliwholecellkb.settings

from ecoliwholecellkb_project.public.models import (Gene, Molecule, Location,
Comment, ProteinMonomers, Rna, Metabolite, ProteinComplex, ProteinComplexModified,
ProteinMonomerModified, RnaModified, RelationStoichiometry,
ProteinComplexReactionRelation,ProteinComplexModifiedReaction,
ProteinComplexModReactionRelation, ProteinComplexModReactionEnzyme,
ProteinMonomerModifiedReaction, ProteinMonomerModReactionEnzyme,
ProteinMonomerModReactionRelation, RnaModifiedReaction, RnaModReactionEnzyme,
RnaModifiedReactionRelation, MetaboliteReaction, MetaboliteReactionEnzyme,
MetaboliteReactionRelation, MetaboliteBiomass, MetaboliteEquivalentEnzyme,
Chromosome, GeneSplices, GeneAbsolutentPosition, EntryPositiveFloatData, GeneType,
Parameter, Constant)

# Import Biopython for sequence handling
import Bio.Seq

# Load units data from Pint
from units.unit_struct_array import UnitStructArray
from units.unit_registration import Q_


class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """


	def __init__(self):

		# Parse data out of database
		self._defineConstants()

		self._loadProducts() # ADDED: for accessing info from other table 
		self._loadComments() # ADDED: for accessing info from other table 
		self._loadCompartments()
		self._loadMetabolites()
		self._loadGenome()
		self._loadGenes()
		self._loadRnas()
		self._loadProteinMonomers() #not dome
		self._createModifiedForms()
		self._loadRelationStoichiometry() # ADDED: for accessing info from other table 
		self._loadComplexes() 
		self._loadReactions()

		## Keep separate
		self._loadBiomassFractions() # Build hacked constants - need to add these to SQL database still
		self._loadConstants()
		self._loadParameters()
		self._loadHacked() 		# Build hacked constants - need to add these to the SQL database still
		self._loadComputeParameters()

		# Create data structures for simulation
		#self._buildCompartments()
		#self._buildBulkMolecules()
		#self._buildBiomass()
		self._buildRnaData()
		self._buildMonomerData()
		self._buildRnaIndexToMonomerMapping()
		self._buildMonomerIndexToRnaMapping()
		self._buildConstants()
		self._buildParameters()
		#self._buildAaCounts()
		#self._buildNtCounts()
		#self._buildRnaExpression()
		#self._buildBiomassFractions()

		# Build dependent calculations
		#self._calculateDependentCompartments()



	def _loadHacked(self):
		# New parameters
		self._parameterData['cellWaterMassFraction'] = Q_(0.7, 'water_g / cell_g')
		self._parameterData['cellDryMassFraction'] = Q_(0.3, 'DCW_g / cell_g')

	def _defineConstants(self):
		from collections import OrderedDict
		self._aaWeights = OrderedDict({
			"A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19,
			"G": 75.07, "H": 155.16, "I": 131.18, "K": 146.19, "L": 131.18,
			"M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20,
			"S": 105.09, "T": 119.12, "U": 168.05, "V": 117.15, "W": 204.23,
			"Y": 181.19
		})

		self._waterWeight = Q_(18.02, 'g / mol')
		self._aaWeightsNoWater = OrderedDict([(key, self._aaWeights[key] - self._waterWeight.magnitude) for key in self._aaWeights])

		# Borrowed from BioPython and modified to be at pH 7.2
		self._ntWeights = OrderedDict({ 
			"A": 345.20,
			"C": 321.18,
			"G": 361.20,
			"U": 322.17,
		})


	def _loadProducts(self):
		# Check database access
		self._checkDatabaseAccess(Molecule)

		# Load products
		all_molecules = Molecule.objects.all()

		self._allProducts 	= dict([(i.id, i.product) for i in all_molecules])

		#ADDED for thisType in loadRelationStoichiometry 
		#updated in RNA, monomer, complex, Metabolite, modifiedForm
		self._allProductType	= dict([(i.product, '') for i in all_molecules])


	def _loadComments(self):		
		# Check database access
		self._checkDatabaseAccess(Comment)

		# Load comments
		all_comments = Comment.objects.all()

		self._allComments = dict([(i.id, i.comment_str) for i in all_comments])


	def _loadCompartments(self):
		# Check database access
		self._checkDatabaseAccess(Location)

		# Load data from Django
		all_locations = Location.objects.all()

		self.compartments = []
		self.compIdToAbbrev = {}
		self._dbLocationId = {} # ADDED: for accessing info from other table 

		# Load data 
		for i in all_locations:			
			c = {"id": i.location_id, "abbrev": i.abbreviation}

			self.compartments.append(c)
			self.compIdToAbbrev[c["id"]] = c["abbrev"]

			self._dbLocationId[i.id] = c["abbrev"]	


	def _loadMetabolites(self):
		self.metabolites = []
		
		#biomass		
		biomass = {}
		all_biomass = MetaboliteBiomass.objects.all()	
		if len(all_biomass) <=0:
			raise Exception, "Database Access Error: Cannot access public_MetaboliteBiomass table"
		
		for i in all_biomass:
			temp = {
					"mmol/gDCW":float(i.biomass_concentration),
					"location": self._dbLocationId[i.biomass_location_fk_id]
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
					'id' : self._allProducts[i.equivalent_enzyme_id_fk_id], 
					'location' : self._dbLocationId[i.location_fk_id]
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
					"id": self._allProducts[i.metabolite_id_id].upper(),
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
					"comments": self._allComments[i.comment_fk_id]
				}
			
			if i.has_biomass:
				m["biomassInfo"] = biomass[i.id]

			if i.has_equivalent_enzyme:
				m["equivEnzIds"] = equ_enz[i.id]		
			
			self.metabolites.append(m)
			self._allProductType[self._allProducts[i.metabolite_id_id]] = 'metabolite' #added


	def _loadBiomassFractions(self):

		doublingTime = [100, 60, 40, 30, 24]

		self._cellDryMassCompositionData = numpy.zeros(len(doublingTime),
			dtype = [('doublingTime',				'f'),
					('proteinMassFraction',			'f'),
					('rnaMassFraction',				'f'),
					('dnaMassFraction',				'f'),
					('lipidMassFraction',			'f'),
					('lpsMassFraction',				'f'),
					('mureinMassFraction',			'f'),
					('glycogenMassFraction',		'f'),
					('solublePoolMassFraction',		'f'),
					('inorganicIonMassFraction',	'f')])

		self._cellDryMassCompositionData['doublingTime'] = doublingTime
		self._cellDryMassCompositionData['proteinMassFraction'] = [0.6756756757,
		0.6046511628, 0.5404157044, 0.5304212168, 0.5202312139]
		self._cellDryMassCompositionData['rnaMassFraction'] = [0.1351351351,
		0.1511627907, 0.1778290993, 0.2059282371, 0.2439306358]
		self._cellDryMassCompositionData['dnaMassFraction'] = [0.0513513514,
		0.0348837209, 0.0260969977, 0.0224648986, 0.0211560694]
		self._cellDryMassCompositionData['lipidMassFraction'] = [0.0655363953,
		0.0995149094, 0.1215552785, 0.1146741575, 0.1020727685]
		self._cellDryMassCompositionData['lpsMassFraction'] = [0.0239595424,
		0.0363817948, 0.0444395642, 0.0419238855, 0.0373169261]
		self._cellDryMassCompositionData['mureinMassFraction'] = [0.0176173106,
		0.0267513197, 0.0326761501, 0.0308263864, 0.0274389163]
		self._cellDryMassCompositionData['glycogenMassFraction'] = [0.0176173106,
		0.0267513197, 0.0326761501, 0.0308263864, 0.0274389163]
		self._cellDryMassCompositionData['solublePoolMassFraction'] = [0.0060603548,
		0.009202454, 0.0112405956, 0.0106042769, 0.0094389872]
		self._cellDryMassCompositionData['inorganicIonMassFraction'] = [0.0070469242,
		0.0107005279, 0.0130704601, 0.0123305546, 0.0109755665]

		## Lipids
		lipidIds = ['pe160[c]','pe160[c]','pe161[c]','pe161[c]','pe181[c]',
		'pe181[c]','pg160[c]','pg160[c]','pg161[c]','pg161[c]','pg181[c]',
		'pg181[c]','clpn160[p]','clpn161[p]','clpn181[p]']
		fracOfLipidMass = [0.0920103159, 0.2365979553, 0.0711465906, 0.1829483758,
		 0.0396601262, 0.1019831816, 0.0443064600, 0.0379769658, 0.0342681284,
		  0.0293726815, 0.0190423107, 0.0163219806, 0.0427979361, 0.0330887203,
		   0.0184782712]
		
		lipidIds = [x[:-3].upper() + x[-3:] for x in lipidIds]
		if abs(sum(fracOfLipidMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellLipidFractionData = numpy.zeros(len(lipidIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'f')])
		self._cellLipidFractionData['metaboliteId'] = lipidIds
		self._cellLipidFractionData['massFraction'] = fracOfLipidMass

		## LPS
		lpsIds = ['colipa[o]']
		fracOfLPSMass = [1.]

		lpsIds = [x[:-3].upper() + x[-3:] for x in lpsIds]
		if abs(sum(fracOfLPSMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellLPSFractionData = numpy.zeros(len(lpsIds), 
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'f')])
		self._cellLPSFractionData['metaboliteId'] = lpsIds
		self._cellLPSFractionData['massFraction'] = fracOfLPSMass

		# Murein
		mureinIds = ['murein4p4p[p]','murein3p3p[p]','murein4px4p[p]',
		'murein3px4p[p]','murein4px4px4p[p]']
		fracOfMureinMass = [0.3959811811, 0.0913460423, 0.397005066,
		0.0423905921, 0.0732771184]

		mureinIds = [x[:-3].upper() + x[-3:] for x in mureinIds]
		if abs(sum(fracOfMureinMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellMureinFractionData = numpy.zeros(len(mureinIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'f')])
		self._cellMureinFractionData['metaboliteId'] = mureinIds
		self._cellMureinFractionData['massFraction'] = fracOfMureinMass

		# Glycogen
		glycogenIds = ['glycogen[c]']
		fracOfGlycogenMass = [1.]

		glycogenIds = [x[:-3].upper() + x[-3:] for x in glycogenIds]
		if abs(sum(fracOfGlycogenMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellGlycogenFractionData = numpy.zeros(len(glycogenIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'f')])
		self._cellGlycogenFractionData['metaboliteId'] = glycogenIds
		self._cellGlycogenFractionData['massFraction'] = fracOfGlycogenMass

		# Soluble pool
		solublePoolIds = ['ptrc[c]','spmd[c]','accoa[c]','coa[c]','succoa[c]',
		'malcoa[c]','nad[c]','nadh[c]','nadp[c]','nadph[c]','fad[c]','thf[c]',
		'mlthf[c]','5mthf[c]','thmpp[c]','q8h2[c]','2dmmql8[c]','mql8[c]',
		'pydx5p[c]','hemeO[c]','pheme[c]','sheme[c]','enter[c]','gthrd[c]',
		'adocbl[c]','udcpdp[c]','10fthf[c]','chor[c]','amet[c]','ribflv[c]']
		fracOfSolublePoolMass = [0.3483925406,0.1161308469,0.0261156573,
		0.0148516941,0.0098435030,0.0030810913,0.1374440284,0.0034413294,
		0.0096012716,0.0288430300,0.0203218325,0.0115004920,0.0118120079,
		0.0118904381,0.0109525963,0.0189109720,0.0182880179,0.0186518206,
		0.0063575867,0.0217069101,0.0179120227,0.0235678621,0.0173654264,
		0.0079446556,0.0409685431,0.0059412634,0.0122269562,0.0058139964,
		0.0103601428,0.0097614647]

		solublePoolIds = [x[:-3].upper() + x[-3:] for x in solublePoolIds]
		if abs(sum(fracOfSolublePoolMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellSolublePoolFractionData = numpy.zeros(len(solublePoolIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'f')])
		self._cellSolublePoolFractionData['metaboliteId'] = solublePoolIds
		self._cellSolublePoolFractionData['massFraction'] = fracOfSolublePoolMass

		# Inorganic ions
		inorganicIonIds = ['k[c]','nh4[c]','mg2[c]','ca2[c]','fe2[c]','fe3[c]',
		'cu2[c]','mn2[c]','mobd[c]','cobalt2[c]','zn2[c]','cl[c]','so4[c]','pi[c]']
		fracInorganicIonMass = [0.6592084209,0.0203462209,0.0180351501,
		0.0180295557,0.0378534452,0.0378534452,0.0191129730,0.0165239120,
		0.0481057645,0.0177255636,0.0192281995,0.0157765944,0.0361161681,
		0.0360845869]

		inorganicIonIds = [x[:-3].upper() + x[-3:] for x in inorganicIonIds]
		if abs(sum(fracInorganicIonMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellInorganicIonFractionData = numpy.zeros(len(inorganicIonIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'f')])
		self._cellInorganicIonFractionData['metaboliteId'] = inorganicIonIds
		self._cellInorganicIonFractionData['massFraction'] = fracInorganicIonMass

	def _loadGenome(self):
		self.translationTable = 11 # E. coli is 11
		
		all_seq = Chromosome.objects.all()
		genome = ''
		for i in all_seq:
			genome = i.sequence
			break
		self.genomeSeq = genome


	def _loadGenes(self):

		self.genes = []
		self.rnas = []
		self.proteins = []
		self._geneDbIds = {} # ADDED: for rnas and monomers

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
			self._geneDbIds[i.id] = i.frame_id # Added for rnas and monomers

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
				g["rnaId"] = self._allProducts[i.productname_id]

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
				"monomerId": None,
				"type": None # Filled in later
			}
			if g["type"] == "mRNA":
				r["name"] = g["name"] + " [RNA]" # else: need to check name in the RNAs file
				r["monomerId"] = self._allProducts[i.productname_id]
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
			self._allProductType[r["id"]] = 'rna' #added

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
				tmp = dict([(x, 0) for x in self._aaWeights])

				for aa in tmp: 
					tmp[aa] = p["seq"].count(aa)
				p["aaCount"] = numpy.array([tmp["A"], tmp["R"], tmp["N"], tmp["D"], tmp["C"],
								tmp["E"], tmp["Q"], tmp["G"], tmp["H"], tmp["I"],
								tmp["L"], tmp["K"], tmp["M"], tmp["F"], tmp["P"],
								tmp["U"], tmp["S"], tmp["T"], tmp["W"], tmp["Y"], tmp["V"]
								])

				water = 18.02
				aaWeights = {}
				for k in self._aaWeights: aaWeights[k] = self._aaWeights[k] - water
				p["mw"] = water
				for aa in p["seq"]: p["mw"] += aaWeights[aa]

				self.proteins.append(p)
				self._allProductType[p["id"]] = 'protein' #added


	def _loadRnas(self):
		#RnaModified
		rnamodified = {}
		all_rnamodified = RnaModified.objects.all()
		if len(all_rnamodified) <=0:
			raise Exception, "Database Access Error: Cannot access public_RnaModified table"
		for i in all_rnamodified:
			if i.unmodified_rna_fk_id not in rnamodified:
				rnamodified[i.unmodified_rna_fk_id] = []
			rnamodified[i.unmodified_rna_fk_id].append(str(self._allProducts[i.rna_mod_id]))	

		#rna
		all_rna = Rna.objects.all()
		if len(all_rna) <=0:
			raise Exception, "Database Access Error: Cannot access public_Rna table"

		# rnaId -> location index in self.rnas
		rnaLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self.rnas)])
		
		for i in all_rna:
			# RNA
			r = {
				"id": self._allProducts[i.frame_id_id],
				"name": i.name,
				"geneId": self._geneDbIds[i.gene_fk_id],
				"location": self._dbLocationId[i.location_fk_id],
				"modifiedForms": [],
				"comments": self._allComments[i.comment_fk_id],
				"type":self.genes[i.gene_fk_id]['type'],
				}
		
			if int(i.is_modified):	
				r["modifiedForms"] = rnamodified[i.id]

			self.rnas[rnaLookup[r["id"]]]["name"] = r["name"]
			self.rnas[rnaLookup[r["id"]]]["location"] = r["location"]
			self.rnas[rnaLookup[r["id"]]]["modifiedForms"] = r["modifiedForms"]
			self.rnas[rnaLookup[r["id"]]]["comments"] = r["comments"]
			self.rnas[rnaLookup[r["id"]]]["type"] = r["type"]


	def _loadProteinMonomers(self):
		
		#ProteinMonomerModified
		monomermod = {}
		all_monomermod = ProteinMonomerModified.objects.all()
		if len(all_monomermod) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinMonomerModified table"
		for i in all_monomermod:
			if i.unmodified_protein_monomer_fk_id not in monomermod:
				monomermod[i.unmodified_protein_monomer_fk_id] = []
			monomermod[i.unmodified_protein_monomer_fk_id].append(str(self._allProducts[i.protein_monomer_mod_id]))

		#ProteinMonomers
		all_monomers = ProteinMonomers.objects.all()
		if len(all_monomers) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinMonomers table"

		# monomerId -> location index in self.proteins
		protLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self.proteins)])

		for i in all_monomers:

			# Monomer
			p = {
				"id": self._allProducts[i.frame_id_id],
				"name": i.name,
				"geneId": self._geneDbIds[i.gene_fk_id],
				"location": self._dbLocationId[i.location_fk_id],
				"modifiedForms": [],
				"comments": self._allComments[i.comment_fk_id]
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

	def _createModifiedForms(self):
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
					self._allProductType[rNew["id"]] = 'rna' #added		
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
					self._allProductType[pNew["id"]] = 'protein' #added		
					proteinsToAppend.append(pNew)

		self.proteins.extend(proteinsToAppend)

	def _loadRelationStoichiometry(self):

		self.allRelationStoichiometry = {}
		
		all_RelationStoichiometry = RelationStoichiometry.objects.all()
		if len(all_RelationStoichiometry) <=0:
			raise Exception, "Database Access Error: Cannot access public_RelationStoichiometry table"

		for i in all_RelationStoichiometry:
			thisType = self._allProductType[self._allProducts[i.reactant_fk_id]]
			self.allRelationStoichiometry[i.id] = { "coeff": float(i.coefficient), "location": self._dbLocationId[i.location_fk_id], "molecule": self._allProducts[i.reactant_fk_id], "form": "mature", "type":  thisType}


	def _loadComplexes(self):
		
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
			complexMod[i.unmodified_protein_complex_fk_id].append(str(self._allProducts[i.protein_complex_mod_id]))	

		#proteinComplexes
		all_complex = ProteinComplex.objects.all()
		if len(all_complex) <=0:
			raise Exception, "Database Access Error: Cannot access public_ProteinComplex table"

		protNew = []		
		for i in all_complex:
			p = {
				"id": self._allProducts[i.protein_complex_id],
				"name": i.name,
				"modifiedForms": [],
				"unmodifiedForm": None,
				"location": self._dbLocationId[i.location_fk_id],
				"composition": [],
				"dir": int(i.reaction_direction),
				"formationProcess": "Complexation",
				"seq": "",
				"aaCount": numpy.zeros(21),
				"ntCount": numpy.zeros(4),
				"mw": -1,
				"geneId": "",
				"rnaId": "",
				"comments": self._allComments[i.comment_fk_id]
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

		self._createModifiedForms()
		
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

	def _loadReactions(self):

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
			enz[i.metabolite_reaction_fk_id].append(str(self._allProducts[i.enzyme_fk_id])) 
		
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

	def _loadConstants(self):
		self._checkDatabaseAccess(Constant)
		all_constant = Constant.objects.all()
		self._constantData = {}
		for c in all_constant:
			self._constantData[c.name] = Q_(c.value, c.units)


	def _loadParameters(self):
		self._checkDatabaseAccess(Parameter)
		all_parameter = Parameter.objects.all()
		self._parameterData = {}
		for p in all_parameter:
			self._parameterData[p.name] = Q_(p.value, p.units)


	def _loadComputeParameters(self):
		self._parameterData['avgCellToInitalCellConvFactor'] = Q_(numpy.exp(numpy.log(2) * self._parameterData['avgCellCellCycleProgress']), 'dimensionless')
		self._parameterData['avgCellDryMassInit'] = self._parameterData['avgCellDryMass'] / self._parameterData['avgCellToInitalCellConvFactor']
		self._parameterData['avgCellWaterMass'] = (self._parameterData['avgCellDryMass'] / self._parameterData['cellDryMassFraction']) * self._parameterData['cellWaterMassFraction']
		self._parameterData['avgCellWaterMassInit'] = self._parameterData['avgCellWaterMass'] / self._parameterData['avgCellToInitalCellConvFactor']

	## -- Build functions -- ##

	def _buildCompartments(self):
		self.compartments 	= self._compartmentData

	def _buildBulkMolecules(self):
		size = len(self._metaboliteData)*len(self._compartmentData) + len(self._rnaData) + len(self._proteinMonomerData)

		self.bulkMolecules = numpy.zeros(size,
			dtype = [("moleculeId", 		"a50"),
					("mass",				"f"),
					("isMetabolite",		"bool"),
					("isRna",				"bool"),
					("isProteinMonomer",	"bool"),
					("isWater",				"bool"),
					("isModifiedForm",		"bool")])

		# Set metabolites
		lastMetaboliteIdx = len(self._metaboliteData) * len(self._compartmentData)
		self.bulkMolecules['moleculeId'][0:lastMetaboliteIdx] = ['{}[{}]'.format(idx,c)
											for c in self._compartmentData['compartmentAbbreviation']
											for idx in self._metaboliteData['id']
											]

		self.bulkMolecules['mass'][0:lastMetaboliteIdx]		= [self._metaboliteData['mw7.2'][i]
											for j in range(len(self._compartmentData))
											for i in range(len(self._metaboliteData['id']))
											]

		self.bulkMolecules['isMetabolite'][0:lastMetaboliteIdx] = [True]*len(self._metaboliteData) * len(self._compartmentData)

		for i,mid in enumerate(self.bulkMolecules['moleculeId']):
			if 'H2O[' == mid[:4]:
				self.bulkMolecules['isWater'][i] = True

		# Set RNA
		lastRnaIdx = len(self._rnaData) + lastMetaboliteIdx
		self.bulkMolecules['moleculeId'][lastMetaboliteIdx:lastRnaIdx] = ['{}[{}]'.format(idx, self._rnaData['location'][i]) for i,idx in enumerate(self._rnaData['id'])]
		self.bulkMolecules['mass'][lastMetaboliteIdx:lastRnaIdx] = self._rnaData['mw']
		self.bulkMolecules['isRna'][lastMetaboliteIdx:lastRnaIdx] = [True]*len(self._rnaData)

		# Set protein monomers
		lastProteinMonomerIdx = len(self._proteinMonomerData) + lastRnaIdx
		self.bulkMolecules['moleculeId'][lastRnaIdx:lastProteinMonomerIdx] = ['{}[{}]'.format(idx, self._proteinMonomerData['location'][i]) for i,idx in enumerate(self._proteinMonomerData['id'])]
		self.bulkMolecules['mass'][lastRnaIdx:lastProteinMonomerIdx] = self._proteinMonomerData['mw']
		self.bulkMolecules['isProteinMonomer'][lastRnaIdx:lastProteinMonomerIdx] = [True]*len(self._proteinMonomerData)

		# Add units to values
		units = {"moleculeId" : None, "mass" : "g / mol", "isMetabolite" : None, "isRna" : None, "isProteinMonomer" : None, "isModifiedForm" : None, 'isWater' : None}
		self.bulkMolecules = UnitStructArray(self.bulkMolecules, units)

	def _buildAaCounts(self): # TODO: ask nick/derek about deleting these methods
		self.proteinMonomerAACounts = self._aaCountData

	def _buildNtCounts(self): # TODO: ask nick/derek about deleting these methods
		self.rnaNTCounts = self._ntCountData

	def _buildRnaExpression(self):
		normalizedRnaExpression = self._rnaData['expression'] / numpy.sum(self._rnaData['expression'])
		self.rnaExpression = Q_(normalizedRnaExpression, 'dimensionless')

	def _buildBiomass(self):
		units = {'metaboliteId' : None,
				'biomassFlux' : 'mmol / (DCW_g*hr)'}
		self.coreBiomass 		= UnitStructArray(self._coreBiomassData, units)
		self.wildtypeBiomass 	= UnitStructArray(self._wildtypeBiomassData, units)

	def _buildBiomassFractions(self):
		units = {
		'doublingTime' : 'min',
		'proteinMassFraction' : None,
		'rnaMassFraction' : None,
		'dnaMassFraction' : None,
		'lipidMassFraction' : None,
		'lpsMassFraction' : None,
		'mureinMassFraction' : None,
		'glycogenMassFraction' : None,
		'solublePoolMassFraction' : None,
		'inorganicIonMassFraction' : None
		}
		self.cellDryMassComposition = UnitStructArray(self._cellDryMassCompositionData, units)
		self.cellLipidFraction = self._cellLipidFractionData
		self.cellLPSFractionData = self._cellLPSFractionData
		self.cellMureinFractionData = self._cellMureinFractionData
		self.cellGlycogenFractionData = self._cellGlycogenFractionData
		self.cellSolublePoolFractionData = self._cellSolublePoolFractionData
		self.cellInorganicIonFractionData = self._cellInorganicIonFractionData

	def _buildRnaData(self):
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location'])
			for rna in self.rnas]

		rnaDegRates = numpy.log(2) / numpy.array(
			[rna['halfLife'] for rna in self.rnas]
			) # TODO: units

		rnaLens = numpy.array([len(rna['seq']) for rna in self.rnas])

		ntCounts = numpy.array([
			(rna['seq'].count('A'), rna['seq'].count('U'),
				rna['seq'].count('C'), rna['seq'].count('G'))
			for rna in self.rnas
			])

		expression = numpy.array([rna['expression'] for rna in self.rnas])

		synthProb = expression * (
			numpy.log(2) / self._parameterData['cellCycleLen'].to('s').magnitude
			+ rnaDegRates)

		synthProb /= synthProb.sum()

		mws = numpy.array([rna['mw'] for rna in self.rnas])

		size = len(rnaIds)

		is23S = numpy.zeros(size, dtype = numpy.bool)
		is16S = numpy.zeros(size, dtype = numpy.bool)
		is5S = numpy.zeros(size, dtype = numpy.bool)

		for rnaIndex, rna in enumerate(self.rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is23S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is16S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is5S[rnaIndex] = True

		# TODO: Add units
		self.rnaData = numpy.zeros(
			size,
			dtype = [
				('id', 'a50'),
				('synthProb', 'f8'),
				('degRate', 'f8'),
				('length', 'i8'),
				('countsAUCG', '4i8'),
				('mw', 'f8'),
				('isMRna', 'bool'),
				('isMiscRna', 'bool'),
				('isRRna', 'bool'),
				('isTRna', 'bool'),
				('isRRna23S', 'bool'),
				('isRRna16S', 'bool'),
				('isRRna5S', 'bool')
				]
			)

		self.rnaData['id'] = rnaIds
		self.rnaData['synthProb'] = synthProb
		self.rnaData['degRate'] = rnaDegRates
		self.rnaData['length'] = rnaLens
		self.rnaData['countsAUCG'] = ntCounts
		self.rnaData['mw'] = mws
		self.rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in self.rnas]
		self.rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in self.rnas]
		self.rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in self.rnas]
		self.rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in self.rnas]
		self.rnaData['isRRna23S'] = is23S
		self.rnaData['isRRna16S'] = is16S
		self.rnaData['isRRna5S'] = is5S


	def _buildMonomerData(self):
		monomers = [protein for protein in self.proteins 
			if protein['unmodifiedForm'] is None
			and len(protein['composition']) == 0]

		ids = ['{}[{}]'.format(protein['id'], protein['location'])
			for protein in monomers]

		rnaIds = []

		for protein in monomers:
			rnaId = protein['rnaId']

			rnaLocation = None
			for rna in self.rnas:
				if rna['id'] == rnaId:
					rnaLocation = rna['location']
					break

			rnaIds.append('{}[{}]'.format(
				rnaId,
				rnaLocation
				))

		lengths = []
		aaCounts = []

		for protein in monomers:
			sequence = protein['seq']

			counts = []

			for aa in self._aaWeights.viewkeys(): # TODO: better way to get AA ids?
				counts.append(
					sequence.count(aa)
					)

			lengths.append(len(sequence))
			aaCounts.append(counts)

		mws = numpy.array([protein['mw'] for protein in monomers])

		size = len(rnaIds)

		nAAs = len(aaCounts[0])

		# TODO: Add units
		self.monomerData = numpy.zeros(
			size,
			dtype = [
				('id', 'a50'),
				('rnaId', 'a50'),
				('length', 'i8'),
				('aaCounts', '{}i8'.format(nAAs)),
				('mw', 'f8'),
				]
			)

		self.monomerData['id'] = ids
		self.monomerData['rnaId'] = rnaIds
		self.monomerData['length'] = lengths
		self.monomerData['aaCounts'] = aaCounts
		self.monomerData['mw'] = mws


	def _buildRnaIndexToMonomerMapping(self):
		self.rnaIndexToMonomerMapping = numpy.array([numpy.where(x == self.rnaData["id"])[0][0] for x in self.monomerData["rnaId"]])


	def _buildMonomerIndexToRnaMapping(self):
		self.monomerIndexToRnaMapping = numpy.array([numpy.where(x == self.monomerData["rnaId"])[0][0] for x in self.rnaData["id"] if len(numpy.where(x == self.monomerData["rnaId"])[0])])


	def _buildConstants(self):
		self.constants = self._constantData
		self.__dict__.update(self.constants)


	def _buildParameters(self):
		self.parameters = self._parameterData
		self.__dict__.update(self.parameters)

## -- Calculate dependent variables -- ##
	def _calculateDependentCompartments(self):
		self.nCompartments 	= len(self.compartments)

## -- Utility functions -- ##
	def _checkDatabaseAccess(self, table):
		if len(table.objects.all()) <= 0:
			raise Exception, "Database Acess Error: Cannot access public_{} table".format(table.__name__.lower())


	def _calculateRnaWeight(self, seq):
		return sum(self._ntWeights[x] for x in seq) - (len(seq) - 1) * 17.01


	def _calculatePeptideWeight(self, seq):
		return sum(self._aaWeightsNoWater[x] for x in seq) + self._waterWeight.to('g/mol').magnitude


	def _calcNucleotideCount(self, seq):
		return numpy.array([seq.count(x) for x in self._ntWeights])


	def _calculateAminoAcidCount(self, seq):
		return numpy.array([seq.count(x) for x in self._aaWeights])
