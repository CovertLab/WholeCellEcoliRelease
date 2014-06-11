#!/usr/bin/env python

"""
KnowledgeBase for Ecoli

Whole-cell knowledge base ecoli

@author: Sajia Akhter
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/04/2014
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/14/2014
"""
from __future__ import division
import numpy # TODO: change to import numpy as np
import collections

import os
import sys
import itertools

# Set Django environmental variable
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb_project.ecoliwholecellkb.settings'

import wholecell.utils.config
sys.path.append(str(os.path.expanduser(wholecell.utils.config.KNOWLEDGEBASE_PACKAGE_DIR)))
import ecoliwholecellkb_project.ecoliwholecellkb.settings

from ecoliwholecellkb_project.public.models import (MoleculeType, Gene, Molecule, Location,
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


AMINO_ACID_1_TO_3_ORDERED = collections.OrderedDict(( # TOKB
	("A", "ALA-L[c]"), ("R", "ARG-L[c]"), ("N", "ASN-L[c]"), ("D", "ASP-L[c]"),
	("C", "CYS-L[c]"), ("E", "GLU-L[c]"), ("Q", "GLN-L[c]"), ("G", "GLY[c]"),
	("H", "HIS-L[c]"), ("I", "ILE-L[c]"), ("L", "LEU-L[c]"), ("K", "LYS-L[c]"),
	("M", "MET-L[c]"), ("F", "PHE-L[c]"), ("P", "PRO-L[c]"), ("S", "SER-L[c]"),
	("T", "THR-L[c]"), ("W", "TRP-L[c]"), ("Y", "TYR-L[c]"), ("U", "SEC-L[c]"),
	("V", "VAL-L[c]")
	))


AMINO_ACID_WEIGHTS = { # TOKB
	"A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19,
	"G": 75.07, "H": 155.16, "I": 131.18, "K": 146.19, "L": 131.18,
	"M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20,
	"S": 105.09, "T": 119.12, "U": 168.05, "V": 117.15, "W": 204.23,
	"Y": 181.19
	}

RNA_TYPE_ORDER = { '23srRNA' : 0, '16srRNA' : 1, '5srRNA' : 2, 'tRNA' : 3, 'mRNA' : 4, 'miscRNA' : 5,'Protein' : 6, 'Metabolite' : 7}

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
		self._loadRelationStoichiometry() #  Need to be called before any reaction loading
		self._loadModificationReactions() # Need to be called before rna/protein/complexes modified forms
		self._loadRnas()
		self._loadProteinMonomers()
		self._loadComplexes() 
		self._loadReactions()

		self._calcMolecularWeightFromRxn()		


		## Keep separate
		self._loadBiomassFractions() # Build hacked constants - need to add these to SQL database still
		self._loadConstants()
		self._loadParameters()
		self._loadHacked() 		# Build hacked constants - need to add these to the SQL database still
		self._loadComputeParameters()
		'''
		# Create data structures for simulation
		self._buildCompartments()
		self._buildBulkMolecules()
		self._buildUniqueMolecules()
		self._buildBiomass()
		self._buildRnaData()
		self._buildMonomerData()
		self._buildRnaIndexToMonomerMapping()
		self._buildMonomerIndexToRnaMapping()
		self._buildConstants()
		self._buildParameters()
		self._buildRnaExpression()
		self._buildBiomassFractions()

		self._buildComplexationMatrix()
		self._buildMetabolism()
		'''
		# Build dependent calculations
		#self._calculateDependentCompartments()

		# delete all auxiliary objects
		self._deleteAuxiliaryVariables()

	def _loadHacked(self):
		# New parameters
		self._parameterData['cellWaterMassFraction'] = Q_(0.7, 'water_g / cell_g')
		self._parameterData['cellDryMassFraction'] = Q_(0.3, 'DCW_g / cell_g')


	def _defineConstants(self):
		self._aaWeights = collections.OrderedDict()

		for singleLetterName in AMINO_ACID_1_TO_3_ORDERED.viewkeys():
			self._aaWeights[singleLetterName] = AMINO_ACID_WEIGHTS[singleLetterName]

		self._waterWeight = Q_(18.02, 'g / mol')
		self._aaWeightsNoWater = collections.OrderedDict([
			(key, self._aaWeights[key] - self._waterWeight.magnitude) for key in self._aaWeights
			])

		# Borrowed from BioPython and modified to be at pH 7.2
		self._ntWeights = collections.OrderedDict({ 
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
		self._checkDatabaseAccess(MoleculeType)
		all_types = MoleculeType.objects.all()
		types  	= dict([(i.id, i.molecule_type) for i in all_types])

		self._allProductType	= dict([(i.product, str(types[i.molecule_type_fk_id])) for i in all_molecules])


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

		self.compartmentList = []
		self._compIdToAbbrev = {}
		self._dbLocationId = {} # ADDED: for accessing info from other table 

		# Load data 
		for i in all_locations:			
			c = {"id": i.location_id, "abbrev": i.abbreviation}

			self.compartmentList.append(c)
			self._compIdToAbbrev[c["id"]] = c["abbrev"]

			self._dbLocationId[i.id] = c["abbrev"]	


	def _loadMetabolites(self):
		self._metabolites = []
		
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
			
			self._metabolites.append(m)
			###self._allProductType[self._allProducts[i.metabolite_id_id]] = 'metabolite' #need to delete

	def _loadRelationStoichiometry(self):

		self._allRelationStoichiometry = {}
		
		all_RelationStoichiometry = RelationStoichiometry.objects.all()
		if len(all_RelationStoichiometry) <=0:
			raise Exception, "Database Access Error: Cannot access public_RelationStoichiometry table"

		for i in all_RelationStoichiometry:
			thisType = self._allProductType[self._allProducts[i.reactant_fk_id]]
			self._allRelationStoichiometry[i.id] = { "coeff": float(i.coefficient), "location": self._dbLocationId[i.location_fk_id], "molecule": self._allProducts[i.reactant_fk_id], "form": "mature", "type":  thisType}

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
		lipidIds = ['pe160[c]','pe160[p]','pe161[c]','pe161[p]','pe181[c]',
		'pe181[p]','pg160[c]','pg160[p]','pg161[c]','pg161[p]','pg181[c]',
		'pg181[p]','clpn160[p]','clpn161[p]','clpn181[p]']
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
		lpsIds = ['colipa[e]']
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

		self._genes = []
		self._geneDbIds = {} # ADDED: for rnas and monomers: will be deleted at the end of DB loading
		
		#genetype
		genetypes = {}
		self._checkDatabaseAccess(GeneType)
		all_genetypes = GeneType.objects.all()
		for i in all_genetypes:
			genetypes[i.id] = i.type_gene
		
		#Gene
		self._checkDatabaseAccess(Gene)		
		all_genes = Gene.objects.all()

		self._expression = {} #addded for RNAs: need to reorganize later
		self._half_life = {} #addded for RNAs: need to reorganize later

		for i in all_genes:
			self._geneDbIds[i.id] = i.frame_id # Added for rnas and monomers
			self._expression[i.frame_id] = i.expression_id #addded for RNAs: need to reorganize later
			self._half_life[i.frame_id] = i.half_life_id #addded for RNAs: need to reorganize later

			g = {
				"id": i.frame_id,
				"name": str(i.name),
				"symbol": i.symbol,
				"type": genetypes[i.typegene_id],
				"coordinate": int(i.coordinate) - 1, # The coordinates we're given are 1 indexed.
				"length": int(i.length),
				"direction": i.direction,
				"seq": "",
				"rnaId": None,
				"monomerId": None
			}
			g["name"] = g["name"].replace("\\","")
			if g["direction"] == "f":
				g["direction"] = '+'
				g["seq"] = self.genomeSeq[(g["coordinate"]): (g["coordinate"] + g["length"])]
			else:
				g["direction"] = '-'
				g["seq"] = Bio.Seq.Seq(self.genomeSeq[(g["coordinate"] - g["length"] + 1): (g["coordinate"] + 1)]).reverse_complement().tostring()

			if g["type"] == "mRNA":
				g["monomerId"] = self._allProducts[i.productname_id]
				g["rnaId"] = g["id"] + "_RNA" ## added for mRNA in loadRNAs()
			else:
				g["rnaId"] = self._allProducts[i.productname_id]

			self._genes.append(g)
	

	def _loadModificationReactions(self):

		self._modificationReactions = []
		self._rnaModReactionDbIds = {} # Added for rnas: will be deleted at the end of DB loading	
		self._proteinModReactionDbIds = {} # Added for monomers: will be deleted at the end of DB loading		
		self._complexModReactionDbIds = {} # Added for complexes: will be deleted at the end of DB loading		

		# modified RNA rxns		
		relation = {}
		self._checkDatabaseAccess(RnaModifiedReactionRelation)		
		all_relation = RnaModifiedReactionRelation.objects.all()
		for i in all_relation:
			if i.rna_mod_reaction_fk_id not in relation:
				relation[i.rna_mod_reaction_fk_id] = []
			relation[i.rna_mod_reaction_fk_id].append(i.reactant_relation_id) 
		
		##
		enz = {}
		self._checkDatabaseAccess(RnaModReactionEnzyme)		
		all_enz = RnaModReactionEnzyme.objects.all()
		for i in all_enz:
			if i.reaction_fk_id not in enz:
				enz[i.reaction_fk_id] = []
			enz[i.reaction_fk_id].append(str(self._allProducts[i.reaction_enzyme_fk_id])) 
		
		##
		self._checkDatabaseAccess(RnaModifiedReaction)		
		all_reaction = RnaModifiedReaction.objects.all()

		for i in all_reaction:

			if i.rna_mod_fk_id not in self._rnaModReactionDbIds:
				self._rnaModReactionDbIds[i.rna_mod_fk_id] = [] 
			self._rnaModReactionDbIds[i.rna_mod_fk_id].append(i.reaction_id)


			r = {
					"id": i.reaction_id,
					"process": "rna",
					"ec": i.ec,
					"dir": int(i.reaction_direction),
					"stoichiometry": [],
					"catBy": []
				}

			if i.id in enz:
				r["catBy"] = enz[i.id]
			else:
				r["catBy"] = None

			if r["ec"] == None: r["ec"] = ""
			
			if i.id not in relation:
				raise Exception, "%s RNA has no reaction" % i.reaction_id
			for temp in relation[i.id]:
				t = self._allRelationStoichiometry[temp]
				#t["molecule"] = t["molecule"].upper() # need to check why .upper()
				r["stoichiometry"].append(t)
	
			self._modificationReactions.append(r)

		# modified monomers rxns		
		relation = {}
		self._checkDatabaseAccess(ProteinMonomerModReactionRelation)		
		all_relation = ProteinMonomerModReactionRelation.objects.all()
		for i in all_relation:
			if i.reaction_fk_id not in relation:
				relation[i.reaction_fk_id] = []
			relation[i.reaction_fk_id].append(i.reactant_relation_id) 
		
		##
		enz = {}
		self._checkDatabaseAccess(ProteinMonomerModReactionEnzyme)		
		all_enz = ProteinMonomerModReactionEnzyme.objects.all()
		for i in all_enz:
			if i.reaction_fk_id not in enz:
				enz[i.reaction_fk_id] = []
			enz[i.reaction_fk_id].append(str(self._allProducts[i.reaction_enzyme_fk_id])) 
		
		##
		self._checkDatabaseAccess(ProteinMonomerModifiedReaction)		
		all_reaction = ProteinMonomerModifiedReaction.objects.all()

		for i in all_reaction:

			if i.protein_monomer_mod_fk_id not in self._proteinModReactionDbIds:
				self._proteinModReactionDbIds[i.protein_monomer_mod_fk_id] = []
			self._proteinModReactionDbIds[i.protein_monomer_mod_fk_id].append(i.reaction_id) 

			r = {
					"id": i.reaction_id,
					"process": "monomer",
					"ec": i.ec,
					"dir": int(i.reaction_direction),
					"stoichiometry": [],
					"catBy": []
				}

			if i.id in enz:
				r["catBy"] = enz[i.id]
			else:
				r["catBy"] = None

			if r["ec"] == None: r["ec"] = ""
			
			if i.id not in relation:
				raise Exception, "%s RNA has no reaction" % i.reaction_id
			for temp in relation[i.id]:
				t = self._allRelationStoichiometry[temp]
				#t["molecule"] = t["molecule"].upper() # need to check why .upper()
				r["stoichiometry"].append(t)
	
			self._modificationReactions.append(r)

		# modified complexes rxns		
		relation = {}
		self._checkDatabaseAccess(ProteinComplexModReactionRelation)		
		all_relation = ProteinComplexModReactionRelation.objects.all()
		for i in all_relation:
			if i.complex_mod_reaction_fk_id not in relation:
				relation[i.complex_mod_reaction_fk_id] = []
			relation[i.complex_mod_reaction_fk_id].append(i.reactant_relation_id) 
		
		##
		enz = {}
		self._checkDatabaseAccess(ProteinComplexModReactionEnzyme)		
		all_enz = ProteinComplexModReactionEnzyme.objects.all()
		for i in all_enz:
			if i.complex_mod_reaction_fk_id not in enz:
				enz[i.complex_mod_reaction_fk_id] = []
			enz[i.complex_mod_reaction_fk_id].append(str(self._allProducts[i.reaction_enzyme_fk_id])) 
		
		##
		self._checkDatabaseAccess(ProteinComplexModifiedReaction)		
		all_reaction = ProteinComplexModifiedReaction.objects.all()

		for i in all_reaction:
			if i.protein_complex_mod_fk_id not in self._complexModReactionDbIds:
				self._complexModReactionDbIds[i.protein_complex_mod_fk_id] = []
			self._complexModReactionDbIds[i.protein_complex_mod_fk_id].append(i.reaction_id) 

			r = {
					"id": i.reaction_id,
					"process": "complex",
					"ec": i.ec,
					"dir": int(i.reaction_direction),
					"stoichiometry": [],
					"catBy": []
				}

			if i.id in enz:
				r["catBy"] = enz[i.id]
			else:
				r["catBy"] = None

			if r["ec"] == None: r["ec"] = ""
			
			if i.id not in relation:
				raise Exception, "%s RNA has no reaction" % i.reaction_id
			for temp in relation[i.id]:
				t = self._allRelationStoichiometry[temp]
				#t["molecule"] = t["molecule"].upper() # need to check why .upper()
				r["stoichiometry"].append(t)
	
			self._modificationReactions.append(r)
		

	def _loadModifiedRnas(self):

		self._modifiedRnas = []
		
		#RnaModified
		rnamodified = {}
		self._checkDatabaseAccess(RnaModified)		
		all_rnamodified = RnaModified.objects.all()
		for i in all_rnamodified:
			rMod = {
				"id": self._allProducts[i.rna_mod_id],
				"name": i.name,
				"location": self._dbLocationId[i.location_fk_id],
				"comments": self._allComments[i.comment_fk_id],
				"reactionId" : self._rnaModReactionDbIds[i.id],
				"mw" : [], 	
				#"unmodifiedForm" : self._allProducts[i.unmodified_rna_fk.frame_id_id] #need to check why gives error!
				"unmodifiedForm" : i.unmodified_rna_fk_id # This is the FK of RNA; Will be updated on _loadRnas()
				}
			
			self._modifiedRnas.append(rMod)
			
			#store for _rnas
			if i.unmodified_rna_fk_id not in rnamodified:
				rnamodified[i.unmodified_rna_fk_id] = []
			rnamodified[i.unmodified_rna_fk_id].append(str(self._allProducts[i.rna_mod_id]))	
			
		return rnamodified

	def _loadRnas(self):

		self._rnas = []
		rnamodified = self._loadModifiedRnas()
		rnaDbIds = {}
		#rna
		self._checkDatabaseAccess(Rna)		
		all_rna = Rna.objects.all()
		
		### geneId -> location index in self._genes
		geneLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self._genes)])
		
		#EntryPositiveFloatData
		posData = {}
		self._checkDatabaseAccess(EntryPositiveFloatData)		
		all_posData = EntryPositiveFloatData.objects.all()
		for i in all_posData:
			posData[i.id] = float(i.value)

		for i in all_rna:
			rnaDbIds[i.id] = self._allProducts[i.frame_id_id]
			gene_frame_id = self._geneDbIds[i.gene_fk_id]

			# RNA
			r = {
				"id": self._allProducts[i.frame_id_id],
				"name": i.name,
				"location": self._dbLocationId[i.location_fk_id],
				"comments": self._allComments[i.comment_fk_id],

				#from other tables
				"modifiedForms": [],	
				"monomerId": None,
				"geneId": gene_frame_id,
				"type": self._genes[geneLookup[gene_frame_id]]["type"],
				"expression": posData[self._expression[gene_frame_id]], #TODO
				"halfLife": posData[self._half_life[gene_frame_id]],	#TODO
								
				#need to calculate				
				"seq": "",
				"ntCount": [],
				"mw": [0, 0, 0, 0, 0, 0, 0, 0]
				}
		
			if int(i.is_modified):	
				if i.id not in rnamodified:
					raise Exception, "%s RNA has no modified form" % i.frame_id_id
				r["modifiedForms"] = rnamodified[i.id]

			if r["type"] == "mRNA":
				r["monomerId"] = self._genes[geneLookup[r["geneId"]]]["monomerId"]
			
			gene_seq = self._genes[geneLookup[r["geneId"]]]["seq"]
			r["seq"] = Bio.Seq.Seq(gene_seq, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()
			r["ntCount"] = numpy.array([r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")])
			weight = 345.20 * r["ntCount"][0] + 321.18 * r["ntCount"][1] + 361.20 * r["ntCount"][2] + 322.17 * r["ntCount"][3] - (len(r["seq"]) - 1) * 17.01
			index = self._whichRna(r['id'], r['type'])
			r["mw"][index] = weight 

						
			self._rnas.append(r)
	
			# TODO from DEREK: Uncomment when Nick has fixed json formatting
			# if type(r["halfLife"]) == dict:
			# 	if r["halfLife"]["units"] != "day":
			# 		raise Exception, "Unknown unit!"
			# 	r["halfLife"] = r["halfLife"]["value"] * 24.0 * 60.0 * 60.0


		##update FK of modified RNAs
		for i in self._modifiedRnas:
			i["unmodifiedForm"] = rnaDbIds[i["unmodifiedForm"]]

		#ADD mRNAs from the GENE table
		for g in self._genes:
			if g["type"] == "mRNA":
				r = {
					"id": g["id"] + "_RNA",
					"name": g["name"] + " [RNA]",
					"location": self._compIdToAbbrev["CCO-CYTOSOL"],
					"comments": "",

					#from other tables
					"modifiedForms": [],		
					"monomerId": g["monomerId"],
					"geneId": g["id"],
					"type": g["type"],
					"expression": posData[self._expression[g["id"]]], #TODO
					"halfLife": posData[self._half_life[g["id"]]],	#TODO
								
					#need to calculate				
					"seq": "",
					"ntCount": [],
					"mw": [0, 0, 0, 0, 0, 0, 0, 0]
				}

				r["seq"] = Bio.Seq.Seq(g["seq"], Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()
				r["ntCount"] = numpy.array([r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")])
				weight = 345.20 * r["ntCount"][0] + 321.18 * r["ntCount"][1] + 361.20 * r["ntCount"][2] + 322.17 * r["ntCount"][3] - (len(r["seq"]) - 1) * 17.01
				index = self._whichRna(r['id'], r['type'])
				r["mw"][index] = weight 
			
				self._rnas.append(r)


	def _whichRna(self, rnaId, rnaType):
		if rnaType == 'miscRNA':
			return RNA_TYPE_ORDER['miscRNA']
		if rnaType == 'tRNA':
			return RNA_TYPE_ORDER['tRNA']
		if rnaType == 'mRNA':
			return RNA_TYPE_ORDER['mRNA']
		if rnaType == "rRNA" and rnaId.startswith("RRL"):
			return RNA_TYPE_ORDER['23srRNA']
		if rnaType == "rRNA" and rnaId.startswith("RRS"):
			return RNA_TYPE_ORDER['16srRNA']
		if rnaType == "rRNA" and rnaId.startswith("RRF"):
			return RNA_TYPE_ORDER['5srRNA']
 
	def _loadModifiedProteinMonomers(self):

		self._modifiedProteins = []
				
		proteinModified = {}
		self._checkDatabaseAccess(ProteinMonomerModified)		
		all_proteinModified = ProteinMonomerModified.objects.all()
		for i in all_proteinModified:
			pMod = {
				"id": self._allProducts[i.protein_monomer_mod_id],
				"name": i.name,
				"location": self._dbLocationId[i.location_fk_id],
				"comments": self._allComments[i.comment_fk_id],
				"reactionId" : None,
				"mw" : [], 	
				#"unmodifiedForm" : self._allProducts[i.unmodified_protein_monomer_fk.frame_id_id] #need to check why gives error!
				"unmodifiedForm" : i.unmodified_protein_monomer_fk_id # This is the FK of RNA; Will be updated on _loadRnas()
				}
			if i.id in self._proteinModReactionDbIds:
				pMod["reactionId"] = self._proteinModReactionDbIds[i.id]
			
			self._modifiedProteins.append(pMod)
			
			#store for _proteins
			if i.unmodified_protein_monomer_fk_id not in proteinModified:
				proteinModified[i.unmodified_protein_monomer_fk_id] = []
			proteinModified[i.unmodified_protein_monomer_fk_id].append(str(self._allProducts[i.protein_monomer_mod_id]))	
			
		return proteinModified

	def _loadProteinMonomers(self):
		
		self._proteins = []
		monomermod = self._loadModifiedProteinMonomers()
		proteinDbIds = {}
		
		### geneId -> location index in self._genes
		geneLookup = dict([(x[1]["id"], x[0]) for x in enumerate(self._genes)])

		#genesplices
		genesplices = {}
		self._checkDatabaseAccess(GeneSplices)		
		all_genesplices = GeneSplices.objects.all()
		for i in all_genesplices:
			genesplices[self._geneDbIds[i.gene_id]] = {
									'start1':int(i.start1),
									'stop1':int(i.stop1), 									
									'start2':int(i.start2),
									'stop2':int(i.stop2)
								}		
		#GeneAbsolutentPosition
		genePos = {}
		self._checkDatabaseAccess(GeneAbsolutentPosition)		
		all_genePos = GeneAbsolutentPosition.objects.all()	
		for i in all_genePos:
			genePos[self._geneDbIds[i.gene_id]] = {'pos':int(i.abs_nt_pos),'old':i.old,'new':i.new}

		#ProteinMonomers
		self._checkDatabaseAccess(ProteinMonomers)		
		all_monomers = ProteinMonomers.objects.all()
		
		for i in all_monomers:
			proteinDbIds[i.id] = self._allProducts[i.frame_id_id]
			gene_frame_id = self._geneDbIds[i.gene_fk_id]
			# Monomer
			p = {
				"id": self._allProducts[i.frame_id_id],
				"name": i.name,
				"geneId": gene_frame_id,
				"location": self._dbLocationId[i.location_fk_id],
				"modifiedForms": [],
				"comments": self._allComments[i.comment_fk_id],
				"seq": "",
				"aaCount": numpy.zeros(21),
				"mw": -1,
				"rnaId": self._genes[geneLookup[gene_frame_id]]["rnaId"]		
			}

			if int(i.is_modified):	#TODO: Check after update monomer_modified by Nick
				if i.id in monomermod: 
					p["modifiedForms"] = monomermod[i.id]
				else:
					raise Exception, "modified Monomer Absent %s" % p["id"]
					#print p["id"], i.id, i.name

				
			
			if gene_frame_id in genesplices:
				baseSequence = Bio.Seq.Seq("", Bio.Alphabet.IUPAC.IUPACUnambiguousDNA())
				baseSequence = baseSequence + self.genomeSeq[genesplices[gene_frame_id]['start1']-1:genesplices[gene_frame_id]['stop1']] + self.genomeSeq[genesplices[gene_frame_id]['start2']-1:genesplices[gene_frame_id]['stop2']]

				if self._genes[geneLookup[gene_frame_id]]["direction"] == "-":
					baseSequence = baseSequence.reverse_complement()
				baseSequence = baseSequence.tostring()

			else:
				baseSequence = self._genes[geneLookup[gene_frame_id]]['seq'] 

			p["seq"] = Bio.Seq.Seq(baseSequence, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).translate(table = self.translationTable).tostring()
			
			if gene_frame_id in genePos:
				pos = genePos[gene_frame_id]['pos']
				before = genePos[gene_frame_id]['old']
				after = genePos[gene_frame_id]['new']
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

			self._proteins.append(p)

		##update FK of modified proteins
		for i in self._modifiedProteins:
			i["unmodifiedForm"] = proteinDbIds[i["unmodifiedForm"]]


	def _loadModifiedProteinComplexes(self):

		self._modifiedComplexes = []
				
		complexModified = {}
		self._checkDatabaseAccess(ProteinComplexModified)		
		all_complexModified = ProteinComplexModified.objects.all()
		for i in all_complexModified:
			pMod = {
				"id": self._allProducts[i.protein_complex_mod_id],
				"name": i.name,
				"location": self._dbLocationId[i.location_fk_id],
				"comments": self._allComments[i.comment_fk_id],
				"reactionId" : None,
				"mw" : [], 	
				#"unmodifiedForm" : self._allProducts[i.unmodified_protein_complex_fk.frame_id_id] #need to check why gives error!
				"unmodifiedForm" : i.unmodified_protein_complex_fk_id # This is the FK of RNA; Will be updated on _loadRnas()
				}
			if i.id in self._complexModReactionDbIds:
				pMod["reactionId"] = self._complexModReactionDbIds[i.id]
			
			self._modifiedComplexes.append(pMod)
			
			#store for _proteinComplexes
			if i.unmodified_protein_complex_fk_id not in complexModified:
				complexModified[i.unmodified_protein_complex_fk_id] = []
			complexModified[i.unmodified_protein_complex_fk_id].append(str(self._allProducts[i.protein_complex_mod_id]))	
			
		return complexModified

	def _loadComplexes(self):
		
		self._complexationReactions = []
		self._proteinComplexes = []
		complexMod = self._loadModifiedProteinComplexes()
		complexDbIds = {}

		##reaction		
		relation = {}
		self._checkDatabaseAccess(ProteinComplexReactionRelation)		
		all_relation = ProteinComplexReactionRelation.objects.all()
		for i in all_relation:
			if i.protein_complex_fk_id not in relation:
				relation[i.protein_complex_fk_id] = []
			relation[i.protein_complex_fk_id].append(i.reactant_relation_id) 

		#proteinComplexes
		self._checkDatabaseAccess(ProteinComplex)		
		all_complex = ProteinComplex.objects.all()

		for i in all_complex:
			complexDbIds[i.id] = self._allProducts[i.protein_complex_id]
			p = {
				"id": self._allProducts[i.protein_complex_id],
				"reactionId": self._allProducts[i.protein_complex_id]+'_RXN',
				"name": i.name,
				"modifiedForms": [],
				"location": self._dbLocationId[i.location_fk_id],
				"mw": [],
				"comments": self._allComments[i.comment_fk_id]
			}
			if i.modified_form:	
				if i.id in complexMod:				
					p["modifiedForms"] = complexMod[i.id]
				else:
					raise Exception, "modification form absent for complex" 
			
			self._proteinComplexes.append(p)

			#reaction
			r = {
				"id": p["reactionId"],
				"process": "complexation",
				"dir": int(i.reaction_direction),
				"stoichiometry": []
				}

			if i.id not in relation:
				raise Exception, "%s RNA has no reaction" % i.protein_complex_id
			for temp in relation[i.id]:
				t = self._allRelationStoichiometry[temp]
				#t["molecule"] = t["molecule"].upper() # need to check why .upper()
				r["stoichiometry"].append(t)
				
			self._complexationReactions.append(r)

		##update FK of modified complexes
		for i in self._modifiedComplexes:
			i["unmodifiedForm"] = complexDbIds[i["unmodifiedForm"]]


	def _loadReactions(self):

		self._reactions = []
		
		##		
		relation = {}
		self._checkDatabaseAccess(MetaboliteReactionRelation)		
		all_relation = MetaboliteReactionRelation.objects.all()
		for i in all_relation:
			if i.metabolite_reaction_fk_id not in relation:
				relation[i.metabolite_reaction_fk_id] = []
			relation[i.metabolite_reaction_fk_id].append(i.reactant_relation_id) 
		
		##
		enz = {}
		self._checkDatabaseAccess(MetaboliteReactionEnzyme)		
		all_enz = MetaboliteReactionEnzyme.objects.all()
		for i in all_enz:
			if i.metabolite_reaction_fk_id not in enz:
				enz[i.metabolite_reaction_fk_id] = []
			enz[i.metabolite_reaction_fk_id].append(str(self._allProducts[i.enzyme_fk_id])) 
		
		##
		self._checkDatabaseAccess(MetaboliteReaction)		
		all_metaboliteReaction = MetaboliteReaction.objects.all()
		

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
				t = self._allRelationStoichiometry[temp]
				t["molecule"] = t["molecule"].upper()
				r["stoichiometry"].append(t)

			self._reactions.append(r)
					

	def _calcMolecularWeightFromRxn(self):
		
		productsPosition = {}
		reactions = {}
		complexReactionLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._complexationReactions)])
		modificationReactionLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._modificationReactions)])
		complexLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._proteinComplexes)])
		rnaModLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._modifiedRnas)])

		# No reaction available for calculating MW
		notComputableList = ['ACECITLY-CPLX','CPLX0-2','ENTB-CPLX','MONOMER0-2863','TAP-GLN',
					'TAR-GLN','TRG-GLN','TSR-GLN','MONOMER0-1842','MONOMER0-1843','MONOMER0-2',
					'MONOMER0-3','MONOMER0-4119','MONOMER0-4195','MONOMER0-4196','MONOMER0-4198',
					'OCTANOYL-ACP','OX-FERREDOXIN','PGLYCEROLTRANSII-MONOMER','PHOSPHASERDECARB-ALPHA-MONOMER',
					'PHOSPHASERDECARB-BETA-MONOMER','PHOSPHO-CPXR','PHOSPHO-DCUR-MONOMER','PHOSPHO-KDPE',
					'PHOSPHO-OMPR-MONOMER','PHOSPHO-RCSB','PHOSPHO-UHPA','SAMDC-ALPHA-MONOMER','SAMDC-BETA-MONOMER',
					#not possible to calculate without RXN of above
					'PHOSPHO-OMPR','CPLX0-7885','ENTMULTI-CPLX','TRG-GLU','PROTEIN-NRIP','CPLX0-7748',
					'PHOSPHASERDECARB-DIMER','CPLX0-7795','SAMDECARB-CPLX','CPLX0-7754','CPLX0-7884','TAP-GLU',
					'CPLX0-263','CPLX0-2901','TSR-GLU','CPLX0-7721','PHOSPHO-BARA-HIS','PHOSPHASERDECARB-CPLX',
					'CPLX0-7978','TRG-GLUME','TAP-GLUME','TSR-GLUME','MONOMER0-1',
					#loop rxn proteinComplex
					'BCCP-CPLX', 'ACETYL-COA-CARBOXYLMULTI-CPLX',
					#loop rxn modproteincomplex
					'PTSI-PHOSPHORYLATED', 'CPLX0-8002', 'CPLX0-8004', 'BCCP-BIOTIN-CO2', 'TAR-GLU', 'TAR-GLUME'
					]
		
		#order: complex, modComplex, modmonomer, modrna 
		index = 0
				
		for i in self._proteinComplexes:
			if i['id'] in notComputableList:	continue
			productsPosition[i['id']] = index
			reactions[i['id']] = [self._complexationReactions[complexReactionLookUp[i['reactionId']]]['stoichiometry']]
			index = index + 1	
		
		for i in self._modifiedComplexes:
			if i['id'] in notComputableList:	continue
			productsPosition[i['id']] = index
			reactions[i['id']] = []
			r = i['reactionId'][len(i['reactionId'])-1] #last reaction; considering only 1 reaction
			#for r in i['reactionId']:
			sch = self._modificationReactions[modificationReactionLookUp[r]]['stoichiometry']
			reactions[i['id']].append(sch)
			index = index + 1
		
		#modMonomer needs to calculate for protein complexes
		monomers_mod = ['AMINOMETHYLDIHYDROLIPOYL-GCVH','LIPOYL-GCVH', 'PHOSPHO-TORS850', 'PHOSPHO-TORASP', 'PHOSPHO-TORS']
		for i in self._modifiedProteins:
			if i['id'] in notComputableList:	continue
			if i['id'] not in monomers_mod: 	continue
			productsPosition[i['id']] = index
			reactions[i['id']] = []
			r = i['reactionId'][len(i['reactionId'])-1] #last reaction; considering only 1 reaction
			#for r in i['reactionId']:
			sch = self._modificationReactions[modificationReactionLookUp[r]]['stoichiometry']
			reactions[i['id']].append(sch)
			index = index + 1	
			
		for i in self._modifiedRnas:
			if i['id'] in notComputableList:	continue
			productsPosition[i['id']] = index
			reactions[i['id']] = []
			r = i['reactionId'][len(i['reactionId'])-1] #last reaction; considering only 1 reaction
			#for r in i['reactionId']:
			sch = self._modificationReactions[modificationReactionLookUp[r]]['stoichiometry']
			reactions[i['id']].append(sch)
			index = index + 1		
			

		# buld matrix for equations
		
		met = dict([(x["id"], x['mw7.2']) for x in self._metabolites])
		rna = dict([(x["id"], sum(x['mw'])) for x in self._rnas])
		monomer = dict([(x["id"], x['mw']) for x in self._proteins])
		
		coefficients = []
		values = []
		totalProducts = len(productsPosition)

		for frame_id in productsPosition:	
			for rxn in reactions[frame_id]: #for a single reaction
				tmpCoef = numpy.zeros(totalProducts)
				tmpvalue = 0
				for i in rxn: #for a molecule in a reaction
					if i['molecule'] in productsPosition: #coef
						pos = productsPosition[i['molecule']]
						tmpCoef[pos] = i['coeff']

					else: #already calculated
						if i['type'] == 'metabolite':
							tmpvalue = tmpvalue + met[i['molecule'].upper()] * i['coeff'] * (-1)
						elif i['type'] == 'rna':
							tmpvalue = tmpvalue + rna[i['molecule']] * i['coeff'] * (-1)
						elif i['type'] == 'proteinmonomers':
							tmpvalue = tmpvalue + monomer[i['molecule']] * i['coeff'] * (-1)
						else:
							#print '\''+frame_id+'\',',i['molecule'], i['type']
							raise Exception, "%s unknown molecule while calculating MW" % i['molecule']
				coefficients.append(tmpCoef)
				values.append(tmpvalue)
		
		coefficients = numpy.array(coefficients)
		values = numpy.array(values)
		mw = numpy.linalg.solve(coefficients, values)		

		#update the MW for proteinComplexes and modifiedRnas
		for i in productsPosition:
			if i in rnaModLookUp:
				weight = mw[productsPosition[i]]
				self._modifiedRnas[rnaModLookUp[i]]['mw'] = [weight]
		
		for i in productsPosition:
			if i in complexLookUp:
				weight = mw[productsPosition[i]]
				self._proteinComplexes[complexLookUp[i]]['mw'] = [weight]

		'''
		for i in self._proteinComplexes:
			if i['mw'] <= 0:
				print i['id'],i['mw']
		'''
						
				
	def _calcKCat(self, enzId, vMax, units):
		if enzId == None or vMax == None:
			return numpy.NaN

		if units == "U/mg":
			prot = next((x for x in self._proteins if x["id"] == enzId), None)
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

	def _deleteAuxiliaryVariables(self):
		del self._allComments
		del self._dbLocationId
		del self._compIdToAbbrev
		del self.compartmentList
		del self._geneDbIds
		del self._expression
		del self._half_life
		del self._rnaModReactionDbIds
		del self._proteinModReactionDbIds
		del self._complexModReactionDbIds



	## -- Build functions -- ##

	def _buildCompartments(self):
		self._compartmentData = numpy.zeros(len(self.compartmentList),
			dtype = [('compartmentId','a20'),('compartmentAbbreviation', 'a1')])

		# Load data into structured array
		self._compartmentData['compartmentId']				= [x['id'] for x in self.compartmentList]
		self._compartmentData['compartmentAbbreviation']	= [x['abbrev'] for x in self.compartmentList]
		self.compartments = self._compartmentData
		self.nCompartments 	= len(self.compartmentList)

	def _buildBulkMolecules(self):
		size = len(self._metabolites)*len(self.compartmentList) + len(self._rnas) + len(self._proteins)
		self.bulkMolecules = numpy.zeros(size,
			dtype = [("moleculeId", 		"a50"),
					('compartment',			"a1"),
					("mass",				"f"),
					("isMetabolite",		"bool"),
					("isRnaMonomer",		"bool"),
					("isProteinMonomer",	"bool"),
					("isWater",				"bool"),
					("isComplex",			"bool"),
					("isModified",			"bool")])

		# Set metabolites
		lastMetaboliteIdx = len(self._metabolites) * len(self.compartmentList)
		self.bulkMolecules['moleculeId'][0:lastMetaboliteIdx] = ['{}[{}]'.format(idx,c)
											for c in [x['abbrev'] for x in self.compartmentList]
											for idx in [x['id'] for x in self._metabolites]
											]

		self.bulkMolecules['mass'][0:lastMetaboliteIdx]		= [self._metabolites[i]['mw7.2']
											for j in range(len(self.compartmentList))
											for i in range(len(self._metabolites))
											]

		self.bulkMolecules['isMetabolite'][0:lastMetaboliteIdx] = [True]*len(self._metabolites) * len(self.compartmentList)

		for i,mid in enumerate(self.bulkMolecules['moleculeId']):
			if mid.startswith('H2O['):
				self.bulkMolecules['isWater'][i] = True
				self.bulkMolecules['isMetabolite'][i] = False

		# Set RNA
		lastRnaIdx = len(self._rnas) + lastMetaboliteIdx
		self.bulkMolecules['moleculeId'][lastMetaboliteIdx:lastRnaIdx] = ['{}[{}]'.format(rna['id'], rna['location']) for rna in self._rnas]
		self.bulkMolecules['mass'][lastMetaboliteIdx:lastRnaIdx] = [x['mw'] for x in self._rnas]
		self.bulkMolecules['isRnaMonomer'][lastMetaboliteIdx:lastRnaIdx] = [False if len(x['composition']) else True for x in self._rnas]
		self.bulkMolecules['isComplex'][lastMetaboliteIdx:lastRnaIdx] = [True if len(x['composition']) else False for x in self._rnas]
		self.bulkMolecules['isModified'][lastMetaboliteIdx:lastRnaIdx] = [True if x['unmodifiedForm'] != None else False for x in self._rnas]

		# Set proteins
		lastProteinMonomerIdx = len(self._proteins) + lastRnaIdx
		self.bulkMolecules['moleculeId'][lastRnaIdx:lastProteinMonomerIdx] = ['{}[{}]'.format(protein['id'],protein['location']) for protein in self._proteins]
		self.bulkMolecules['mass'][lastRnaIdx:lastProteinMonomerIdx] = [x['mw'] for x in self._proteins]
		self.bulkMolecules['isModified'][lastRnaIdx:lastProteinMonomerIdx] = [True if x['unmodifiedForm'] != None else False for x in self._proteins]
		self.bulkMolecules['isProteinMonomer'][lastRnaIdx:lastProteinMonomerIdx] = [False if len(x['composition']) else True for x in self._proteins]
		self.bulkMolecules['isComplex'][lastRnaIdx:lastProteinMonomerIdx] = [True if len(x['composition']) else False for x in self._proteins]

		# Add units to values
		units = {"moleculeId"	:	None,
			"mass"				:	"g / mol",
			'compartment'		:	None,
			"isMetabolite"		:	None,
			"isRnaMonomer" 		:	None,
			"isProteinMonomer"	:	None,
			"isModified"		:	None,
			'isWater'			:	None,
			'isComplex'			:	None}
		self.bulkMolecules = UnitStructArray(self.bulkMolecules, units)


	def _buildUniqueMolecules(self):
		# TODO: ask Nick about the best way to use the unit struct arrays here
		G_PER_MOL_TO_FG_PER_MOLECULE = 1e15 / 6.022e23

		self.uniqueMoleculeDefinitions = {
			'activeRnaPoly' : {
				'rnaIndex' : 'i8',
				'requiredACGU' : '4i8',
				'assignedACGU' : '4i8',
				},
			'activeRibosome' : {
				'proteinIndex' : 'i8',
				'requiredAAs' : '20i8',
				'assignedAAs' : '20i8'
				},
			'activeDnaPolymerase' : {
				'chromosomeLocation' : 'i64',
				'directionIsPositive' : 'bool'
				}
			}

		rnaPolyComplexSubunits = [
			"EG10893-MONOMER", "EG10893-MONOMER",			# 2 sub-units are required
			"RPOB-MONOMER", "RPOC-MONOMER", "RPOD-MONOMER"
			]

		rnaPolyComplexMass = sum(
			protein['mw'] for protein in self._proteins
			if protein['id'] in rnaPolyComplexSubunits
			)

		# TODO: This is a bad hack that works because in the fitter
		# I have forced expression to be these subunits only
		ribosomeSubunits = [
			"RRLA-RRNA", "RRSA-RRNA", "RRFA-RRNA"
			]

		ribosomeMass = sum(
			rna['mw'] for rna in self._rnas
			if rna['id'] in ribosomeSubunits
			)

		self.uniqueMoleculeMasses = numpy.zeros(
			shape = len(self.uniqueMoleculeDefinitions),
			dtype = [
				('moleculeId', 'a50'),
				('massMetabolite', numpy.float),
				('massRna', numpy.float),
				('massProtein', numpy.float),
				]
			)

		self.uniqueMoleculeMasses[0] = (
			'activeRnaPoly',
			0,
			0,
			rnaPolyComplexMass * G_PER_MOL_TO_FG_PER_MOLECULE
			)
		self.uniqueMoleculeMasses[1] = (
			'activeRibosome',
			0,
			ribosomeMass * G_PER_MOL_TO_FG_PER_MOLECULE,
			0
			)
		self.uniqueMoleculeMasses[2] = (
			'activeDnaPolymerase',
			0,
			0,
			0
			)

		# TODO: units
		# TODO: make this logic better overall


	def _buildRnaExpression(self):
		normalizedRnaExpression = numpy.zeros(sum(1 for x in self._rnas if x['unmodifiedForm'] == None),
			dtype = [('rnaId',		'a50'),
					('expression',	'f'),
					('isMRna',		'bool'),
					('isMiscRna',	'bool'),
					('isRRna',		'bool'),
					('isTRna',		'bool'),
					('isRRna23S',	'bool'),
					('isRRna16S',	'bool'),
					('isRRna5S',	'bool')])

		normalizedRnaExpression['rnaId'] 		= ['{}[{}]'.format(x['id'], x['location']) for x in self._rnas if x['unmodifiedForm'] == None]
		normalizedRnaExpression['expression']	= [x['expression'] for x in self._rnas if x['unmodifiedForm'] == None]
		normalizedRnaExpression['expression']	= normalizedRnaExpression['expression'] / numpy.sum(normalizedRnaExpression['expression'])
		normalizedRnaExpression['isMRna'] = [rna["type"] == "mRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		normalizedRnaExpression['isMiscRna'] = [rna["type"] == "miscRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		normalizedRnaExpression['isRRna'] = [rna["type"] == "rRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		normalizedRnaExpression['isTRna'] = [rna["type"] == "tRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]

		self.rnaExpression = UnitStructArray(normalizedRnaExpression,
			{
			'rnaId'		:	None,
			'expression':	'dimensionless',
			'isMRna'	:	None,
			'isMiscRna'	:	None,
			'isRRna'	:	None,
			'isTRna'	:	None,
			'isRRna23S'	:	None,
			'isRRna16S'	:	None,
			'isRRna5S'	:	None
			})


	def _buildBiomass(self):
		self._coreBiomassData = numpy.zeros(sum(len(x['biomassInfo']['core']) for x in self._metabolites if len(x['biomassInfo']['core'])),
			dtype = [('metaboliteId', 'a50'),
					('biomassFlux', 	'f')])

		self._wildtypeBiomassData = numpy.zeros(sum(len(x['biomassInfo']['wildtype']) for x in self._metabolites if len(x['biomassInfo']['wildtype'])),
			dtype = [('metaboliteId', 'a50'),
					('biomassFlux',		'f')])

		self._coreBiomassData['metaboliteId']	= [
		'{}[{}]'.format(x['id'], x['biomassInfo']['core'][i]['location'])
		for x in self._metabolites if len(x['biomassInfo']['core'])
		for i in range(len(x['biomassInfo']['core']))
		]
		
		self._coreBiomassData['biomassFlux']	= [
		x['biomassInfo']['core'][i]['mmol/gDCW']
		for x in self._metabolites if len(x['biomassInfo']['core'])
		for i in range(len(x['biomassInfo']['core']))
		]
		
		self._wildtypeBiomassData['metaboliteId']	= [
		'{}[{}]'.format(x['id'], x['biomassInfo']['wildtype'][i]['location'])
		for x in self._metabolites if len(x['biomassInfo']['wildtype'])
		for i in range(len(x['biomassInfo']['wildtype']))
		]

		self._wildtypeBiomassData['biomassFlux']	= [
		x['biomassInfo']['wildtype'][i]['mmol/gDCW']
		for x in self._metabolites if len(x['biomassInfo']['wildtype'])
		for i in range(len(x['biomassInfo']['wildtype']))
		]

		units = {'metaboliteId' : None,
				'biomassFlux' : 'mmol / (DCW_g)'}
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
		self.cellLipidFractionData = self._cellLipidFractionData
		self.cellLPSFractionData = self._cellLPSFractionData
		self.cellMureinFractionData = self._cellMureinFractionData
		self.cellGlycogenFractionData = self._cellGlycogenFractionData
		self.cellSolublePoolFractionData = self._cellSolublePoolFractionData
		self.cellInorganicIonFractionData = self._cellInorganicIonFractionData

	def _buildRnaData(self):
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location'])
			for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]

		rnaDegRates = numpy.log(2) / numpy.array(
			[rna['halfLife'] for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
			) # TODO: units

		rnaLens = numpy.array([len(rna['seq']) for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0])

		ntCounts = numpy.array([
			(rna['seq'].count('A'), rna['seq'].count('C'),
				rna['seq'].count('G'), rna['seq'].count('U'))
			for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0
			])

		expression = numpy.array([rna['expression'] for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0])

		synthProb = expression * (
			numpy.log(2) / self._parameterData['cellCycleLen'].to('s').magnitude
			+ rnaDegRates)

		synthProb /= synthProb.sum()

		mws = numpy.array([rna['mw'] for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0])

		size = len(rnaIds)

		is23S = numpy.zeros(size, dtype = numpy.bool)
		is16S = numpy.zeros(size, dtype = numpy.bool)
		is5S = numpy.zeros(size, dtype = numpy.bool)

		for rnaIndex, rna in enumerate(self._rnas):
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
				# TODO: add expression to this table
				('synthProb', 'f8'),
				('degRate', 'f8'),
				('length', 'i8'),
				('countsACGU', '4i8'),
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
		self.rnaData['countsACGU'] = ntCounts
		self.rnaData['mw'] = mws
		self.rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		self.rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		self.rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		self.rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in self._rnas
			if rna['unmodifiedForm'] is None
			and len(rna['composition']) == 0]
		self.rnaData['isRRna23S'] = is23S
		self.rnaData['isRRna16S'] = is16S
		self.rnaData['isRRna5S'] = is5S

		units = {
		'id'		:	None,
		'synthProb' :	'dimensionless',
		'degRate'	:	'1 / s',
		'length'	:	'nucleotide',
		'countsACGU':	'nucleotide',
		'mw'		:	'g / mol',
		'isMRna'	:	None,
		'isMiscRna'	:	None,
		'isRRna'	:	None,
		'isTRna'	:	None,
		'isRRna23S'	:	None,
		'isRRna16S'	:	None,
		'isRRna5S'	:	None
		}


		self.rnaData = UnitStructArray(self.rnaData, units)

	def _buildMonomerData(self):
		monomers = [protein for protein in self._proteins 
			if protein['unmodifiedForm'] is None
			and len(protein['composition']) == 0]

		ids = ['{}[{}]'.format(protein['id'], protein['location'])
			for protein in monomers]

		rnaIds = []

		for protein in monomers:
			rnaId = protein['rnaId']

			rnaLocation = None
			for rna in self._rnas:
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

		self.aaIDs = AMINO_ACID_1_TO_3_ORDERED.values()


		units = {
		'id'		:	None,
		'rnaId'		:	None,
		'length'	:	'amino_acid',
		'aaCounts'	:	'amino_acid',
		'mw'		:	'g / mol'
		}

		self.monomerData = UnitStructArray(self.monomerData, units)


	def _buildRnaIndexToMonomerMapping(self):
		self.rnaIndexToMonomerMapping = numpy.array([numpy.where(x == self.rnaData["id"])[0][0] for x in self.monomerData["rnaId"]])


	def _buildMonomerIndexToRnaMapping(self):
		self.monomerIndexToRnaMapping = numpy.array([numpy.where(x == self.monomerData["rnaId"])[0][0] for x in self.rnaData["id"] if len(numpy.where(x == self.monomerData["rnaId"])[0])])


	def _buildComplexationMatrix(self):
		# Builds a matrix that maps complexes (on the columns) to the correct
		# stoichiometric ratios of subunits (on the rows)

		# Note that this is NOT a valid S matrix for a complexation process; 
		# it is only a graph

		complexesToSubunits = collections.defaultdict(dict)
		subunits = set()

		complexes = complexesToSubunits.viewkeys()

		for molecule in self._proteins + self._rnas:
			composition = molecule['composition']

			if composition:
				complexName = '{}[{}]'.format(molecule['id'], molecule['location'])

				assert complexName not in complexes, 'Duplicate complex ID'

				for subunit in composition:
					coeff = subunit['coeff']

					if coeff > 0: # entry is actually the complex itself
						# Make sure the data makes sense
						assert molecule['id'] == subunit['molecule']
						assert molecule['location'] == subunit['location']
						assert coeff == 1


					else: # entry is a true subunit
						assert coeff % 1 == 0, 'Noninteger subunit stoichiometry'

						subunitName = '{}[{}]'.format(subunit['molecule'], subunit['location'])

						assert subunitName not in complexesToSubunits[complexName], 'Duplicate subunit ID'

						complexesToSubunits[complexName][subunitName] = -coeff

						subunits.add(subunitName)

		complexNames = sorted(complexes)
		subunitNames = sorted(subunits)

		nComplexes = len(complexNames)
		nSubunits = len(subunitNames)

		complexNameToIndex = {complexName:i for i, complexName in enumerate(complexNames)}
		subunitNameToIndex = {subunitName:i for i, subunitName in enumerate(subunitNames)}

		matrix = numpy.zeros((nSubunits, nComplexes), numpy.int64)

		for complexName, subunits in complexesToSubunits.viewitems():
			complexIndex = complexNameToIndex[complexName]

			for subunitName, count in subunits.viewitems():
				subunitIndex = subunitNameToIndex[subunitName]

				matrix[subunitIndex, complexIndex] = count

		self.complexationMatrix = matrix
		self.complexationMatrixComplexIds = complexNames
		self.complexationMatrixSubunitIds = subunitNames


	def _buildMetabolism(self):
		# Build the matrices/vectors for metabolism (FBA)

		# These may be modified/extended later, but should provide the basic
		# data structures

		# Collect reaction information

		allEnzymes = []
		allReversibility = []
		allReactionStoich = []

		molecules = set()

		for reaction in self._reactions:
			assert reaction["process"] == "Metabolism"

			enzymes = reaction['catBy']

			reversible = (reaction['dir'] == 0)

			reactionStoich = {
				'{}[{}]'.format(reactant['molecule'], reactant['location']) : reactant['coeff']
				for reactant in reaction['stoichiometry']
				} 

			allEnzymes.append(enzymes)
			allReversibility.append(reversible)
			allReactionStoich.append(reactionStoich)

			molecules |= reactionStoich.viewkeys()

		nEdges = len(allEnzymes)
		nNodes = len(molecules)

		stoichMatrix = numpy.zeros((nNodes, nEdges))

		# TODO: actually track/annotate enzymes, k_cats

		self.metabolismMoleculeNames = sorted(molecules)

		moleculeNameToIndex = {
			molecule:i
			for i, molecule in enumerate(self.metabolismMoleculeNames)
			}

		# Fill in the basic reaction matrix (as described by Feist)

		for reactionIndex, reactionStoich in enumerate(allReactionStoich):
			for molecule, stoich in reactionStoich.viewitems():
				moleculeIndex = moleculeNameToIndex[molecule]

				stoichMatrix[moleculeIndex, reactionIndex] = stoich

		# Collect exchange reactions

		## First, find anything that looks like an exchange reaction

		exchangeIndexes = numpy.where((stoichMatrix != 0).sum(0) == 1)[0]

		exchangeNames = [
			self.metabolismMoleculeNames[
				numpy.where(stoichMatrix[:, reactionIndex])[0][0]
				]
			for reactionIndex in exchangeIndexes
			]

		## Separate intercellular (sink) vs. extracellular (media) exchange fluxes

		sinkIndexes = []
		sinkNames = []

		externalIndexes = []
		externalNames = []

		for index, name in itertools.izip(exchangeIndexes, exchangeNames):
			if name.endswith('[e]'):
				externalIndexes.append(index)
				externalNames.append(name)

			else:
				sinkIndexes.append(index)
				sinkNames.append(name)

		self.metabolismSinkExchangeReactionIndexes = numpy.array(sinkIndexes)
		self.metabolismSinkExchangeReactionNames = sinkNames

		self.metabolismMediaExchangeReactionIndexes = numpy.array(externalIndexes)
		self.metabolismMediaExchangeReactionNames = externalNames

		# Extend stoich matrix to include intercellular exchange fluxes

		# There are a few possible approaches that need to be discussed/tested
		# Exhaustive: exchange fluxes for every non-extracellular metabolite
		# Simplified: exchange fluxes for every imported extracellular metabolite
		# Rational: actually look at the network structure
		# Empirical: run simulation and see what is used

		# For now, I'm taking the exhaustive approach

		internalNames = [
			moleculeName
			for moleculeName in self.metabolismMoleculeNames
			if not moleculeName.endswith('[e]') and moleculeName not in sinkNames
			]

		internalIndexes = []

		self.metabolismStoichMatrix = numpy.hstack([
			stoichMatrix,
			numpy.zeros((nNodes, len(internalNames)))
			])

		# TODO: decide whether these exchange reactions should be reversible
		allReversibility.extend([True]*len(internalNames))

		for i, name in enumerate(internalNames):
			reactionIndex = nEdges + i

			moleculeIndex = moleculeNameToIndex[name]

			self.metabolismStoichMatrix[moleculeIndex, reactionIndex] = +1

			internalIndexes.append(reactionIndex)

		self.metabolismReversibleReactions = numpy.array(allReversibility, numpy.bool)

		self.metabolismInternalExchangeReactionIndexes = numpy.array(internalIndexes)
		self.metabolismInternalExchangeReactionNames = internalNames


	def _buildConstants(self):
		self.constants = self._constantData
		self.__dict__.update(self.constants)


	def _buildParameters(self):
		self.parameters = self._parameterData
		self.parameters["fracActiveRibosomes"] = Q_(1.0, "dimensionless")
		self.__dict__.update(self.parameters)

## -- Utility functions -- ##
	def _checkDatabaseAccess(self, table):
		if len(table.objects.all()) <= 0:
			raise Exception, "Database Access Error: Cannot access public_{} table".format(table.__name__.lower())


	def _calculateRnaWeight(self, seq):
		# Starting with NTP molecular weights, subtracting OH for each bond pair
		return sum(self._ntWeights[x] for x in seq) - (len(seq) - 1) * 17.01


	def _calculatePeptideWeight(self, seq):
		return sum(self._aaWeightsNoWater[x] for x in seq) + self._waterWeight.to('g/mol').magnitude


	def _calcNucleotideCount(self, seq):
		return numpy.array([seq.count(x) for x in self._ntWeights])


	def _calculateAminoAcidCount(self, seq):
		return numpy.array([seq.count(x) for x in self._aaWeights])
