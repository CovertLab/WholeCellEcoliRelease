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
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb_project.ecoliwholecellkb.settings'
import ecoliwholecellkb_project.ecoliwholecellkb.settings

import Bio.Seq

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

# Load units data from Pint
from ecoliwholecellkb_project.units.unit_struct_array import UnitStructArray
from ecoliwholecellkb_project.units.unit_registration import Q_


class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """


	def __init__(self):

		# Parse data out of database
		self._defineConstants()
		self._loadProducts() 
		self._loadComments()
		self._loadCompartments()
		self._loadMetabolites()
		self._loadBiomassFractions()
		self._loadGenome()
		self._loadGenes()
		self._loadRnas()
		self._loadProteinMonomers()
		# self.createModifiedForms()
		# self.loadRelationStoichiometry() # ADDED: for accessing info from other table 
		# self.loadComplexes() 
		# self._loadRelationStoichiometry()
		# self._loadReactions()
		self._loadConstants()
		self._loadParameters()

		# Build hacked constants - need to add these to the SQL database still
		self._loadHacked()

		self._loadComputeParameters()



		# Create data structures for simulation
		self._buildCompartments()
		self._buildBulkMolecules()
		self._buildBiomass()
		self._buildRnaData()
		self._buildMonomerData()
		self._buildRnaIndexToMonomerMapping()
		self._buildConstants()
		self._buildParameters()
		self._buildAaCounts()
		self._buildNtCounts()
		self._buildRnaExpression()
		self._buildBiomassFractions()

		# Build dependent calculations
		self._calculateDependentCompartments()



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

		# Initalize structurd array
		self._compartmentData = numpy.zeros(len(all_locations),
			dtype = [('compartmentId','a20'),('compartmentAbbreviation', 'a1')])

		# Load data into structured array
		self._compartmentData['compartmentId']				= [x.location_id for x in all_locations]
		self._compartmentData['compartmentAbbreviation']	= [x.abbreviation for x in all_locations]

		# Save foreign keys for other data parsing
		self._dbLocationId = dict([(x.id, x.abbreviation) for x in all_locations])


	def _loadMetabolites(self):
		## Load Metabolites
		##################

		# Check database access
		self._checkDatabaseAccess(Metabolite)

		# Load data from Django
		all_metabolites = Metabolite.objects.all()

		# TODO: Add compartments to ID strings?
		self._metaboliteData = numpy.zeros(len(all_metabolites),
			dtype = [("id", 		"a50"),
					("name", 		"a200"),
					("feistFormula","a20"),
					("formula7.2",	"a20"),
					("charge7.2",	"i"),
					("mw7.2",		"f"),
					("mediaConc",	"f"),
					("maxExchange",	"f"),
					("fakeMet",		"b"),
					("comments",	"a100")])

		self._metaboliteData['id']				= [self._allProducts[x.metabolite_id_id].upper()
												 for x in all_metabolites]
		self._metaboliteData['name']			= [x.name for x in all_metabolites]
		self._metaboliteData["feistFormula"]	= [x.feist_formula for x in all_metabolites]
		self._metaboliteData["formula7.2"]		= [x.ph_formula for x in all_metabolites]
		self._metaboliteData["charge7.2"]		= [int(x.ph_charge) for x in all_metabolites]
		self._metaboliteData["mw7.2"]			= [float(x.ph_weight) for x in all_metabolites]
		self._metaboliteData["mediaConc"]		= [float(x.media_concentration) for x in all_metabolites]
		self._metaboliteData["maxExchange"]		= [float(x.maximum_exchange_rate)
												for x in all_metabolites]
		self._metaboliteData["fakeMet"]			= [x.fake_metabolite for x in all_metabolites]
		self._metaboliteData["comments"]		= [self._allComments[x.comment_fk_id]
												for x in all_metabolites]

		metabolite_ids_to_string = dict([(x.id,self._allProducts[x.metabolite_id_id].upper())
											 for x in all_metabolites])

		for i in all_metabolites:
			self._allProductType[self._allProducts[i.metabolite_id_id]] = 'metabolite'

		## Load biomass data
		####################

		# Check database access
		self._checkDatabaseAccess(MetaboliteBiomass)

		# Load data from Django
		all_biomass = MetaboliteBiomass.objects.all()

		self._coreBiomassData = numpy.zeros(sum(1 for x in all_biomass if x.is_core),
			dtype = [('metaboliteId', 'a50'),
					('biomassFlux', 	'f'),
					('maintenanceFlux', 'f')])

		self._wildtypeBiomassData = numpy.zeros(sum(1 for x in all_biomass if x.is_wildtype),
			dtype = [('metaboliteId', 'a50'),
					('biomassFlux',		'f'),
					('maintenanceFlux', 'f')])

		self._coreBiomassData['metaboliteId']	= [metabolite_ids_to_string[x.metabolite_id_fk_id]
												+ '[' + self._dbLocationId[x.biomass_location_fk_id] + ']'
												for x in all_biomass if x.is_core]
		self._coreBiomassData['biomassFlux']	= [x.biomass_concentration for x in all_biomass if x.is_core]

		self._wildtypeBiomassData['metaboliteId']	= [metabolite_ids_to_string[x.metabolite_id_fk_id]
													+ '[' + self._dbLocationId[x.biomass_location_fk_id] + ']'
													for x in all_biomass if x.is_wildtype]
		self._wildtypeBiomassData['biomassFlux']	= [x.biomass_concentration for x in all_biomass if x.is_wildtype]

		# TODO: Why is wild type biomass the same length as core biomass? Checked original
		# knowledgebase biomass, they are the same there too!

		## Load equivalent enzyme ID data for pseudo metabolites
		########################################################

		# Check database access
		self._checkDatabaseAccess(MetaboliteEquivalentEnzyme)

		# Load data from Django
		all_equ_enz = MetaboliteEquivalentEnzyme.objects.all()

		self._pseudoMetaboliteEnzymeData = numpy.zeros(len(all_equ_enz),
			dtype = [('metaboliteId', 'a50'),
					('enzyme_id', 'a50')])

		self._pseudoMetaboliteEnzymeData['metaboliteId'] = [metabolite_ids_to_string[x.metabolite_id_fk_id] for x in all_equ_enz]
		self._pseudoMetaboliteEnzymeData['enzyme_id'] = [self._allProducts[x.equivalent_enzyme_id_fk_id] + '[' + self._dbLocationId[x.location_fk_id] + ']' for x in all_equ_enz]


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

		# Load gene types (i.e. mRNA, rRNA, etc.)
		self._checkDatabaseAccess(GeneType)
		all_genetypes = GeneType.objects.all()
		genetypes = dict([(i.id, i.type_gene) for i in all_genetypes])

		# Load gene expression and half life data table
		posData = {}
		all_posData = EntryPositiveFloatData.objects.all()
		if len(all_posData) <=0:
			raise Exception, "Database Access Error: Cannot access public_EntryPositiveFloatData table"
		for i in all_posData:
			posData[i.id] = float(i.value)

		# Load data for genes
		self._checkDatabaseAccess(Gene)
		all_genes = Gene.objects.all()

		self._geneData = numpy.zeros(len(all_genes),
			dtype = [("id", 				"a50"),
					("name", 				"a130"),
					("symbol", 				"a10"),
					("type", 				"a10"),
					("coordinate", 			"i"),
					("length", 				"i"),
					("direction", 			"a1"),
					("sequence", 			"a7110"),
					("rnaId",				"a15"),
					("spliced",				"bool"),
					("aaSubstitution",		"bool")])

		self._geneData["id"]		= [i.frame_id for i in all_genes]
		self._geneData["name"]		= [str(i.name).replace("\\","") for i in all_genes]
		self._geneData["symbol"]	= [i.symbol for i in all_genes]
		self._geneData["type"]		= [genetypes[i.typegene_id] for i in all_genes]
		self._geneData["coordinate"]= [int(i.coordinate) - 1 for i in all_genes]
		self._geneData["length"]	= [int(i.length) for i in all_genes]
		self._geneData["direction"]	= ['+' if i.direction=='f' else '-' for i in all_genes]
		self._geneData["sequence"]	= [self.genomeSeq[(int(i.coordinate) - 1): (int(i.coordinate) - 1 + int(i.length))]
									if i.direction=='f'
									else Bio.Seq.Seq(self.genomeSeq[(int(i.coordinate) - 1 - int(i.length) + 1): (int(i.coordinate) - 1 + 1)]).reverse_complement().tostring()
									for i in all_genes]
		self._geneData["rnaId"] = [i.frame_id if genetypes[i.typegene_id] == "mRNA" else self._allProducts[i.productname_id] for i in all_genes]
		self._geneData["spliced"] = [True if int(i.splices) else False for i in all_genes]
		self._geneData["aaSubstitution"] = [True if int(i.absolute_nt_position) else False for i in all_genes]

		# Data used in creating protein monomers
		self._geneExpressionData = [posData[i.expression_id] for i in all_genes]
		self._geneHalfLifeData = [posData[i.half_life_id] for i in all_genes]
		self._geneMonomerData = [self._allProducts[i.productname_id] if genetypes[i.typegene_id] == 'mRNA' else None for i in all_genes]
		
		# Data used to access gene ids by foreign key
		self._geneDbIds = dict([(i.id, i.frame_id) for i in all_genes])

		# Data used to access gene ids by frame id
		self._geneFrameIdToDbId = dict([(i.frame_id, i.id) for i in all_genes])

		# Save sequence data for translated product. This is non-ideal that this
		# is stored here but its too inconvientent to put it in loadProteinMonomers
		#######################################################################
		# Load gene splices data from frame shifts in translated form
		self._checkDatabaseAccess(GeneSplices)
		all_genesplices = GeneSplices.objects.all()
		genesplices = dict([(i.gene_id,
			{'start1':int(i.start1),'stop1':int(i.stop1),'start2':int(i.start2),
			'stop2':int(i.stop2)}) for i in all_genesplices])	

		# Load substitutions for amino acids in translated sequence for gene product
		self._checkDatabaseAccess(GeneAbsolutentPosition)
		all_genePos = GeneAbsolutentPosition.objects.all()
		genePos = dict([(i.gene_id, {'pos':int(i.abs_nt_pos),
			'old':i.old,'new':i.new}) for i in all_genePos])

		# Generate translated sequences
		self._geneSequenceWithSplicesAndSubstitutions = {}
		for j,seq in enumerate(self._geneData['sequence']):
			geneId = self._geneFrameIdToDbId[self._geneData['id'][j]]
			if self._geneData['spliced'][j]:				
				baseSequence = Bio.Seq.Seq("", Bio.Alphabet.IUPAC.IUPACUnambiguousDNA())
				baseSequence = baseSequence + self.genomeSeq[genesplices[geneId]['start1']-1:genesplices[geneId]['stop1']] + self.genomeSeq[genesplices[geneId]['start2']-1:genesplices[geneId]['stop2']]

				if self._geneData['direction'][j] == "-":
					baseSequence = baseSequence.reverse_complement()
				baseSequence = baseSequence.tostring()

			else:
				baseSequence = seq

			final_sequence = Bio.Seq.Seq(baseSequence, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).translate(table = self.translationTable).tostring()

			if self._geneData['aaSubstitution'][j]:
				pos = genePos[geneId]['pos']
				before = genePos[geneId]['old']
				after = genePos[geneId]['new']
				seqList = list(final_sequence)

				if seqList[pos - 1] != before:
					raise Exception, "Amino acid substitution appears to be incorrect."
				else:
					seqList[pos - 1] = after
				final_sequence = "".join(seqList)

			# Remove * for termation codon
			final_sequence = final_sequence[:final_sequence.find('*')]

			# Save for later useage in loading protein monomers
			self._geneSequenceWithSplicesAndSubstitutions[self._geneData['id'][j]] = final_sequence


	def _loadRnas(self):
		## Create non-coding RNAs #############################################
		#######################################################################
		# self._checkDatabaseAccess(Rna)
		# all_ncrna = Rna.objects.all()

		# # NOTE: This does not load modified forms of these RNAs!
		# self._nonCodingRnaData = numpy.zeros(sum(1 for i in all_ncrna if not int(i.is_modified)),
		# 	dtype = [("id",			"a50"),
		# 			("name",		"a130"),
		# 			("expression",	"float"),
		# 			("halfLife",	"float"),
		# 			("location",	"a1"),
		# 			("sequence",	"a8000"), #TODO: Check
		# 			("mw", 			"f"),
		# 			("geneId", 		"a15"),
		# 			("monomerId",	"a15"),
		# 			("type",		"a10")])
		# self._nonCodingRnaData['id'] = [self._allProducts[i.frame_id_id] for i in all_ncrna if not int(i.is_modified)]
		# self._nonCodingRnaData['name'] = [i.name for i in all_ncrna if not int(i.is_modified)]
		# self._nonCodingRnaData['expression']
		# self._nonCodingRnaData['halfLife']
		# self._nonCodingRnaData['location'] = [self._dbLocationId[i.location_fk_id] for i in all_ncrna if not int(i.is_modified)]
		# self._nonCodingRnaData['sequence']
		# self._nonCodingRnaData['mw'] # Dependent on sequence
		# self._nonCodingRnaData['geneId']
		# self._nonCodingRnaData['monomerId']
		# self._nonCodingRnaData['type']

		## Create mRNAs from genes ############################################
		#######################################################################

		# Load expression data for every gene
		self._checkDatabaseAccess(EntryPositiveFloatData)
		all_posData = EntryPositiveFloatData.objects.all()
		posData = dict([(i.id, i.value) for i in all_posData])

		self._rnaData = numpy.zeros(len(self._geneData),
			dtype = [("id",			"a50"),
					("name",		"a130"),
					("expression",	"float"),
					("halfLife",	"float"),
					("location",	"a1"),
					("sequence",	"a8000"), #TODO: Check
					("mw", 			"f"),
					("geneId", 		"a15"),
					("monomerId",	"a15"),
					("type",		"a10")])

		self._rnaData["id"] = [x for x in self._geneData["rnaId"]]
		self._rnaData["name"] = [x + " [RNA]" for x in self._geneData["name"]]
		self._rnaData["expression"] = self._geneExpressionData
		self._rnaData["halfLife"] = self._geneHalfLifeData
		self._rnaData["location"] = ["c"]*len(self._rnaData)
		self._rnaData["sequence"] = [Bio.Seq.Seq(x, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring() for x in self._geneData["sequence"]]
		self._rnaData["mw"] = [self._calculateRnaWeight(x) for x in self._rnaData["sequence"]]
		self._rnaData["geneId"] = [x for x in self._geneData["id"]]
		self._rnaData["monomerId"] = self._geneMonomerData
		self._rnaData["type"] = self._geneData["type"]

		self._ntCountData = numpy.zeros((len(self._geneData), 4))
		for i in range(len(self._rnaData)):
			self._ntCountData[i,:] = self._calcNucleotideCount(self._rnaData["sequence"][i])

		# TODO: Add
		# "modifiedForms": [],
		# "unmodifiedForm": None,
		# "composition": []

		## Create modified RNAs ###############################################
		#######################################################################

		# TODO: Load modified forms


	def _loadProteinMonomers(self):
		## Load protein monomers ##############################################
		#######################################################################

		# Load data from database
		self._checkDatabaseAccess(ProteinMonomers)
		all_monomers = ProteinMonomers.objects.all()

		self._proteinMonomerData = numpy.zeros(sum(self._geneData["type"] == "mRNA"),
			dtype = [("id",		"a50"),
					("name",	"a130"),
					("location","a1"),
					("sequence","a2400"),
					("mw",		"f"),
					("geneId",	"a50"),
					("rnaId",	"a50")])

		self._proteinMonomerData["id"] = [self._allProducts[i.frame_id_id] for i in all_monomers]
		self._proteinMonomerData["name"] = [i.name for i in all_monomers]
		self._proteinMonomerData["location"] = [self._dbLocationId[i.location_fk_id] for i in all_monomers]
		self._proteinMonomerData["geneId"] = [self._geneDbIds[i.gene_fk_id] for i in all_monomers]
		self._proteinMonomerData["rnaId"] = [self._geneDbIds[i.gene_fk_id] for i in all_monomers]
		self._proteinMonomerData["sequence"] = [self._geneSequenceWithSplicesAndSubstitutions[self._geneDbIds[i.gene_fk_id]] for i in all_monomers]
		self._proteinMonomerData["mw"] = [self._calculatePeptideWeight(x) for x in self._proteinMonomerData['sequence']]

		self._aaCountData = numpy.zeros((sum(self._geneData["type"] == "mRNA"),21))
		for i in range(len(self._proteinMonomerData)):
			self._aaCountData[i,:] = self._calculateAminoAcidCount(self._proteinMonomerData["sequence"][i])

		## Load modified protein monomers #####################################
		#######################################################################

		# TODO: Fill in
		# "modifiedForms": [],
		# "unmodifiedForm": None,
		# "composition": [],
		# "formationProcess": "",


	def _loadRelationStoichiometry(self):

		self.allRelationStoichiometry = {}
		self._checkDatabaseAccess(RelationStoichiometry)
		all_RelationStoichiometry = RelationStoichiometry.objects.all()

		for i in all_RelationStoichiometry:
			thisType = self._allProductType[self._allProducts[i.reactant_fk_id]]
			self.allRelationStoichiometry[i.id] = { "coeff": float(i.coefficient), "location": self._dbLocationId[i.location_fk_id], "molecule": self._allProducts[i.reactant_fk_id], "form": "mature", "type":  thisType}

	def _loadReactions(self):

		self.reactions = []
		self.reactionsExchange = []
		
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

			if r["name"] == None:
				r["name"] = ""
			if r["ec"] == None:
				r["ec"] = ""
			
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

	def _buildAaCounts(self):
		self.proteinMonomerAACounts = self._aaCountData

	def _buildNtCounts(self):
		self.rnaNTCounts = self._ntCountData

	def _buildRnaExpression(self):
		normalizedRnaExpression = self._rnaData['expression'] / numpy.sum(self._rnaData['expression'])
		self.rnaExpression = Q_(normalizedRnaExpression, 'dimensionless')

	def _buildBiomass(self):
		units = {'metaboliteId' : None,
				'biomassFlux' : 'mmol / (DCW_g*hr)',
				'maintenanceFlux' : 'mmol / (DCW_g*hr)'}
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
		rnaIds = ['{}[{}]'.format(id_, location) for id_, location 
			in zip(self._rnaData['id'], self._rnaData['location'])]

		rnaDegRates = numpy.log(2) / self._rnaData['halfLife']

		rnaLens = numpy.array([len(s) for s in self._rnaData['sequence']])

		ntCounts = numpy.array([
			(s.count('A'), s.count('U'), s.count('C'), s.count('G'))
			for s in self._rnaData['sequence']
			])

		synthProb = self._rnaData['expression'] * (
			numpy.log(2) / self._parameterData['cellCycleLen'].to('s').magnitude
			+ numpy.log(2) / self._rnaData['halfLife'])

		synthProb /= synthProb.sum()

		size = len(rnaIds)

		is23S = numpy.zeros(size, dtype = numpy.bool)
		is16S = numpy.zeros(size, dtype = numpy.bool)
		is5S = numpy.zeros(size, dtype = numpy.bool)

		for idx, rna in enumerate(self._rnaData):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is23S[idx] = True
			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is16S[idx] = True
			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is5S[idx] = True

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
		self.rnaData['mw'] = self._rnaData["mw"]
		self.rnaData['isMRna'] = self._rnaData["type"] == "mRNA"
		self.rnaData['isMiscRna'] = self._rnaData["type"] == "miscRNA"
		self.rnaData['isRRna'] = self._rnaData["type"] == "rRNA"
		self.rnaData['isTRna'] = self._rnaData["type"] == "tRNA"
		self.rnaData['isRRna23S'] = is23S
		self.rnaData['isRRna16S'] = is16S
		self.rnaData['isRRna5S'] = is5S


	def _buildMonomerData(self):
		ids = ['{}[{}]'.format(id_, location) for id_, location 
			in zip(self._proteinMonomerData['id'], self._proteinMonomerData['location'])]

		rnaIds = []

		for rnaId in self._proteinMonomerData['rnaId']:
			index = numpy.where(rnaId == self._rnaData['id'])[0][0]

			rnaIds.append('{}[{}]'.format(
				rnaId,
				self._rnaData['location'][index]
				))

		lengths = []
		aaCounts = []

		for sequence in self._proteinMonomerData['sequence']:
			counts = []

			for aa in self._aaWeights.viewkeys(): # TODO: better way to get AA ids?
				counts.append(
					sequence.count(aa)
					)

			lengths.append(len(sequence))
			aaCounts.append(counts)

		size = len(rnaIds)

		nAAs = len(aaCounts[0])

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
		self.monomerData['mw'] = self._proteinMonomerData["mw"]


	def _buildRnaIndexToMonomerMapping(self):
		self.rnaIndexToMonomerMapping = numpy.array([numpy.where(x == self.rnaData["id"])[0][0] for x in self.monomerData["rnaId"]])


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
