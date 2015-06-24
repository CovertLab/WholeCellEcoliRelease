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

import collections
import os
import sys
import itertools
import re

# Import Biopython for sequence handling
import Bio
import Bio.Seq
import warnings
warnings.simplefilter("ignore", Bio.BiopythonWarning)

import numpy as np
import scipy.constants

# Set Django environmental variable
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb_project.ecoliwholecellkb.settings'

import wholecell.utils.config
sys.path.append(str(os.path.expanduser(wholecell.utils.config.KNOWLEDGEBASE_PACKAGE_DIR)))
import ecoliwholecellkb_project.ecoliwholecellkb.settings

from ecoliwholecellkb_project.public.models import *


# Load units data
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils import units

# NOTE: most hard coded constants have been moved to another .py file
#		to keep this file length managable.
from hard_coded_data import *

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	submassNameToIndex = MOLECULAR_WEIGHT_ORDER

	def __init__(self, deleteLoadingData = True):

		# Cache attribute names so we can delete those made by the _load* methods
		defaultAttrs = set(dir(self))

		# Parse data out of database
		self._defineConstants()

		self._loadProducts() # ADDED: for accessing info from other table
		self._loadComments() # ADDED: for accessing info from other table
		self._loadCompartments()
		self._loadMetabolites()
		self._loadPolymerized()
		self._loadGenome()
		self._loadGenes()
		self._loadRelationStoichiometry() #  Need to be called before any reaction loading
		self._loadModificationReactions() # Need to be called before rna/protein/complexes modified forms
		self._loadRnas()
		self._loadProteinMonomers()
		self._loadComplexes()
		self._loadReactions()
		self._loadMetaboliteConcentrations()

		self._calcMolecularWeightFromRxn()

		## Keep separate
		self._loadBiomassFractions() # Build hacked constants - need to add these to SQL database still
		self._loadConstants()
		self._loadParameters()
		self._loadHacked() 		# Build hacked constants - need to add these to the SQL database still
		#self._calcMolecularWeightFromRxn() # Have to call again to calculate MWs of hacked complexes
		self._loadComputeParameters()

		loadedAttrs = set(dir(self)) - defaultAttrs

		self._loadPromoters()
		self._loadTerminators()
		self._loadTranscriptionUnits()
		#self._countATinPromoters()

		# Create data structures for simulation
		self._buildAllMasses() # called early because useful for other builders # Mass # DONE
		self._buildMoleculeGroups() # called early because useful for other builders # State (?) # DONE

		self._buildSequence() # Replication # DONE
		self._buildCompartments() # State (?) # DONE
		self._buildBulkMolecules() # State (?) # DONE
		self._buildBulkChromosome() # State (?) # DONE
		self._buildGeneData() # Replication # DONE
		self._buildRibosomeData() # Translation # SKIP
		self._buildUniqueMolecules() # State (?) # DONE
		self._buildBiomass() # Metabolism (?) # DONE
		self._buildRnaData() # Transcription # DONE
		self._buildMonomerData() # Translation # DONE
		self._buildRnaIndexToMonomerMapping() # Translation (?) # DONE
		self._buildMonomerIndexToRnaMapping() # Translation (?) # DONE
		self._buildRnaIndexToGeneMapping() # Transcription (?) # DONE
		self._buildConstants() # Constants # DONE
		self._buildParameters() # Constants # DONE
		self._buildRnaExpression() # Transcription (*actually goes to raw data)# DONE
		self._buildBiomassFractions() # Metabolism (?)# TODO: move to new KB
		self._buildTranscription() # Transcription # DONE
		self._buildTranslation() # Translation # DONE
		self._buildMetabolitePools() # Metabolism (?)# DONE
		self._buildTrnaData() # Translation (?) # SKIP

		from .complexation import Complexation
		self.complexation = Complexation(self) # DONE

		from .metabolism import Metabolism
		self.metabolism = Metabolism(self) # DONE

		# Build dependent calculations
		#self._calculateDependentCompartments()

		if deleteLoadingData:
			for attr in loadedAttrs:
				delattr(self, attr)


	def _loadHacked(self):
		# New parameters
		self._parameterData['cellWaterMassFraction'] = 0.7
		self._parameterData['cellDryMassFraction'] = 0.3
		self._parameterData['dnaPolymeraseElongationRate'] = 750*units.nt/units.s
		self._parameterData['oriCCenter'] = 3923882*units.nt
		self._parameterData['terCCenter'] = 1607192*units.nt
		self._parameterData['gtpPerTranslation'] = 4.2 # TODO: find a real number
		self._parameterData["fractionChargedTrna"] = 0.8


		# Assumed reaction for producing L-selenocysteine without a tRNA
		# [c]: SER-L + SELNP --> SEC-L + PI + (2)H
		stoich = [{'coeff': -1.0,
				'form': 'mature',
				'location': u'c',
				'molecule': u'SER-L',
				'type': 'metabolite'},
				{'coeff': -1.0,
				'form': 'mature',
				'location': u'c',
				'molecule': u'SELNP',
				'type': 'metabolite'},
				{'coeff': 1.0,
				'form': 'mature',
				'location': u'c',
				'molecule': u'SEC-L',
				'type': 'metabolite'},
				{'coeff': 1.0,
				'form': 'mature',
				'location': u'c',
				'molecule': u'PI',
				'type': 'metabolite'},
				{'coeff': 2.0,
				'form': 'mature',
				'location': u'c',
				'molecule': u'H',
				'type': 'metabolite'}]

		reaction = {
				"id": 'SEC_RXN_HACKED',
				"name": 'L-selenocysteine production reaction (HACKED)',
				"process": "Metabolism",
				"ec": '',
				"dir": 1,
				"stoichiometry": stoich,
				"catBy": [],
				"ub": 1000.,
				"lb": -1000.,
				"kcat": None
			}

		self._reactions.append(reaction)

		# SELNP media exchange
		stoich = [
			{'coeff': 1.0,
			'form': 'mature',
			'location': u'e',
			'molecule': u'SELNP',
			'type': 'metabolite'},
			]

		reaction = {
			"id": 'SELNP_MEDIA_EXCHANGE_HACKED',
			"name": 'Selenium media exchange (HACKED)',
			"process": "Metabolism",
			"ec": '',
			"dir": 0, # presumably reversible
			"stoichiometry": stoich,
			"catBy": [],
			"ub": 1000.,
			"lb": -1000.,
			"kcat": None
			}

		self._reactions.append(reaction)

		# TODO: e->p->c instead of directly to the cytoplasm

		# SELNP import reaction
		stoich = [
			{'coeff': -1.0,
			'form': 'mature',
			'location': u'e',
			'molecule': u'SELNP',
			'type': 'metabolite'},
			{'coeff': 1.0,
			'form': 'mature',
			'location': u'c',
			'molecule': u'SELNP',
			'type': 'metabolite'}
			]

		reaction = {
			"id": 'SELNP_IMPORT_HACKED',
			"name": 'Selenium import reaction (HACKED)',
			"process": "Metabolism",
			"ec": '',
			"dir": 0, # presumably reversible
			"stoichiometry": stoich,
			"catBy": [],
			"ub": 1000.,
			"lb": -1000.,
			"kcat": None
			}

		self._reactions.append(reaction)

		# Fixing localizations for two ribosomal proteins and the complexes they are in
		for prot in self._proteins:
			if prot['id'] == 'EG10877-MONOMER':
				prot['location'] = u'c'
			if prot['id'] == 'EG10876-MONOMER':
				prot['location'] = u'c'

		for comp in self._proteinComplexes:
			if comp['id'] == 'CPLX0-3962':
				comp['location'] = u'c'

		for rxn in self._complexationReactions:
			for molecule in rxn['stoichiometry']:
				if molecule['molecule'] == 'CPLX0-3962':
					molecule['location'] = u'c'
				elif molecule['molecule'] == 'EG10877-MONOMER':
					molecule['location'] = u'c'
				elif molecule['molecule'] == 'EG10876-MONOMER':
					molecule['location'] = u'c'

		# Changing ids of 30S and 50S ribosomal complexes
		for comp in self._proteinComplexes:
			if comp['id'] == 'CPLX0-3953':
				comp['id'] = 'CPLX-30SA'
			elif comp['id'] == 'CPLX0-3962':
				comp['id'] = 'CPLX-50SA'

		for rxn in self._complexationReactions:
			for molecule in rxn['stoichiometry']:
				if molecule['molecule'] == 'CPLX0-3953':
					molecule['molecule'] = 'CPLX-30SA'
					molecule['name'] = '30S ribosomal subunit rrnA'
				elif molecule['molecule'] == 'CPLX0-3962':
					molecule['molecule'] = 'CPLX-50SA'
					molecule['name'] = '50S ribosomal subunit rrnA'

		for rxn in self._complexationReactions:
			if rxn['id'] == 'CPLX0-3953_RXN':
				rxn['id'] = 'CPLX-30SA_RXN'
			elif rxn['id'] == 'CPLX0-3962_RXN':
				rxn['id'] = 'CPLX-50SA_RXN'

		# Add other rrn operons and their formation reactions
		# Ignoring extra 5S rRNA
		remaining16SrRNA = S30_16S_RRNAS[1:]
		letters = [x[3] for x in remaining16SrRNA]
		for idx,rRNA in enumerate(remaining16SrRNA):
			newComplex = {
				'comments': u'',
				'id': u'CPLX-30S{}'.format(letters[idx]),
				'location': u'c',
				'mw': np.zeros(len(MOLECULAR_WEIGHT_ORDER)),
				#'name': u'30S ribosomal subunit rrn{}'.format(letters[idx]),
				'reactionId': u'CPLX-30S{}_RXN'.format(letters[idx])}

			self._proteinComplexes.append(newComplex)

			newStoichiometry = []
			for protein_idx,protein in enumerate(S30_PROTEINS):
				newSubunit = {
					'coeff': -1.*S30_PROTEINS_STOICHIOMETRY[protein_idx],
					'form': 'mature',
					'location': u'c',
					'molecule': protein[:-3],
					'type': 'proteinmonomers'
					}
				newStoichiometry.append(newSubunit)

			newStoichiometry.append({
					'coeff': -1.*S30_16S_RRNAS_STOICHIOMETRY[idx],
					'form': 'mature',
					'location': u'c',
					'molecule': rRNA[:-3],
					'type': 'rna'
					})

			newStoichiometry.append({
					'coeff': 1.,
					'form': 'mature',
					'location': u'c',
					'molecule': 'CPLX-30S{}'.format(letters[idx]),
					'type': 'proteincomplex'
					})

			newComplexationReaction = {
				'dir' : 1,
				'id': u'CPLX-30S{}_RXN'.format(letters[idx]),
				'process' : 'complexation',
				'stoichiometry' : newStoichiometry
				}

			self._complexationReactions.append(newComplexationReaction)

		remaining5SrRNA = S50_5S_RRNAS[1:]
		remaining23SrRNA = S50_23S_RRNAS[1:]
		letters = [x[3] for x in remaining23SrRNA]
		for idx in range(len(remaining23SrRNA)):
			newComplex = {
				'comments': u'',
				'id': u'CPLX-50S{}'.format(letters[idx]),
				'location': u'c',
				'mw': np.zeros(len(MOLECULAR_WEIGHT_ORDER)),
				#'name': u'50S ribosomal subunit rrn{}'.format(letters[idx]),
				'reactionId': u'CPLX-50S{}_RXN'.format(letters[idx])}

			self._proteinComplexes.append(newComplex)

			newStoichiometry = []
			for protein_idx,protein in enumerate(S50_PROTEINS):
				newSubunit = {
					'coeff': -1.*S50_PROTEINS_STOICHIOMETRY[protein_idx],
					'form': 'mature',
					'location': u'c',
					'molecule': protein[:-3],
					'type': 'proteinmonomers'
					}
				newStoichiometry.append(newSubunit)


			for cplx_idx,cplx in enumerate(S50_PROTEIN_COMPLEXES):
				newSubunit = {
					'coeff': -1.*S50_PROTEIN_COMPLEXES_STOICHIOMETRY[cplx_idx],
					'form': 'mature',
					'location': u'c',
					'molecule': cplx[:-3],
					'type': 'proteincomplex'
					}
				newStoichiometry.append(newSubunit)

			newStoichiometry.append({
					'coeff': -1.*S50_23S_RRNAS_STOICHIOMETRY[idx],
					'form': 'mature',
					'location': u'c',
					'molecule': remaining23SrRNA[idx][:-3],
					'type': 'rna'
					})

			newStoichiometry.append({
					'coeff': -1.*S50_5S_RRNAS_STOICHIOMETRY[idx],
					'form': 'mature',
					'location': u'c',
					'molecule': remaining5SrRNA[idx][:-3],
					'type': 'rna'
					})

			newStoichiometry.append({
					'coeff': 1.,
					'form': 'mature',
					'location': u'c',
					'molecule': 'CPLX-50S{}'.format(letters[idx]),
					'type': 'proteincomplex'
					})

			newComplexationReaction = {
				'dir' : 1,
				'id': u'CPLX-30S{}_RXN'.format(letters[idx]),
				'process' : 'complexation',
				'stoichiometry' : newStoichiometry
				}

			self._complexationReactions.append(newComplexationReaction)

	def _defineConstants(self):
		self._aaWeights = collections.OrderedDict()

		for singleLetterName in AMINO_ACID_1_TO_3_ORDERED.viewkeys():
			self._aaWeights[singleLetterName] = None # placeholder

		self._polypeptideEndWeight = None

		self._ntWeights = collections.OrderedDict([
			("A", None),
			("C", None),
			("G", None),
			("U", None),
			])

		self._dntWeights = collections.OrderedDict([
			("A", None),
			("C", None),
			("G", None),
			("T", None),
			])


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

		self._compartmentList = []
		self._compIdToAbbrev = {}
		self._dbLocationId = {} # ADDED: for accessing info from other table

		# Load data
		for i in all_locations:
			c = {"id": i.location_id, "abbrev": i.abbreviation}

			self._compartmentList.append(c)
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


	def _loadPolymerized(self):
		self._polymerized = []

		# Load AAs

		aminoAcidIDtoSingleLetter = {
			value.replace("[c]", ""):key
			for key, value in AMINO_ACID_1_TO_3_ORDERED.viewitems()
			}

		aminoAcidIDtoSortedIndex = {
			value.replace("[c]", ""):index
			for index, value in enumerate(AMINO_ACID_1_TO_3_ORDERED.viewvalues())
			}

		self._polymerizedAA_IDs = [None] * len(aminoAcidIDtoSortedIndex)

		proteinMassKey = MOLECULAR_WEIGHT_KEYS.index("protein")

		for monomer in POLYMERIZED_AMINO_ACID_WEIGHTS:
			entry = {
				"id":monomer["frame id"],
				"mw":monomer["mw"],
				"mass key":proteinMassKey
				}

			self._polymerized.append(entry)

			singleLetter = aminoAcidIDtoSingleLetter[monomer["base molecule"]]

			self._aaWeights[singleLetter] = monomer["mw"]

			index = aminoAcidIDtoSortedIndex[monomer["base molecule"]]

			self._polymerizedAA_IDs[index] = monomer["frame id"]

		# Load peptide end weight (= 1 water)

		entry = {
			"id":POLYPEPTIDE_END_WEIGHT["frame id"],
			"mw":POLYPEPTIDE_END_WEIGHT["mw"],
			"mass key":proteinMassKey,
			}

		self._polymerized.append(entry)

		self._polypeptideEndWeight = POLYPEPTIDE_END_WEIGHT["mw"] # TODO: rename this attribute

		# Load nucleotides

		ntpIDtoSingleLetter = {
			"ATP":"A",
			"CTP":"C",
			"GTP":"G",
			"UTP":"U",
			}

		ntpIDtoSortedIndex = {
			"ATP":0,
			"CTP":1,
			"GTP":2,
			"UTP":3,
			}

		self._polymerizedNT_IDs = [None] * len(ntpIDtoSortedIndex)

		rnaMassKey = MOLECULAR_WEIGHT_KEYS.index("RNA")

		for monomer in POLYMERIZED_NUCLEOTIDE_WEIGHTS:
			entry = {
				"id":monomer["frame id"],
				"mw":monomer["mw"],
				"mass key":rnaMassKey
				}

			self._polymerized.append(entry)

			singleLetter = ntpIDtoSingleLetter[monomer["base molecule"]]

			self._ntWeights[singleLetter] = monomer["mw"]

			index = ntpIDtoSortedIndex[monomer["base molecule"]]

			self._polymerizedNT_IDs[index] = monomer["frame id"]

		# Load RNA end weight (= 1 PPi)

		entry = {
			"id":RNA_END_WEIGHT["frame id"],
			"mw":RNA_END_WEIGHT["mw"],
			"mass key":rnaMassKey,
			}

		self._polymerized.append(entry)

		self._rnaEndWeight = RNA_END_WEIGHT["mw"] # TODO: rename this attribute

		# Load deoxynucleotides

		dntpIDtoSingleLetter = {
			"DATP":"A",
			"DCTP":"C",
			"DGTP":"G",
			"DTTP":"T",
			}

		dntpIDtoSortedIndex = {
			"DATP":0,
			"DCTP":1,
			"DGTP":2,
			"DTTP":3,
			}

		self._polymerizedDNT_IDs = [None] * len(dntpIDtoSortedIndex)

		dnaMassKey = MOLECULAR_WEIGHT_KEYS.index("DNA")

		for monomer in POLYMERIZED_DEOXY_NUCLEOTIDE_WEIGHTS:
			entry = {
				"id":monomer["frame id"],
				"mw":monomer["mw"],
				"mass key":dnaMassKey
				}

			self._polymerized.append(entry)

			singleLetter = dntpIDtoSingleLetter[monomer["base molecule"]]

			self._dntWeights[singleLetter] = monomer["mw"]

			index = dntpIDtoSortedIndex[monomer["base molecule"]]

			self._polymerizedDNT_IDs[index] = monomer["frame id"]


	def _loadRelationStoichiometry(self):

		self._allRelationStoichiometry = {}

		all_RelationStoichiometry = RelationStoichiometry.objects.all()
		if len(all_RelationStoichiometry) <=0:
			raise Exception, "Database Access Error: Cannot access public_RelationStoichiometry table"

		for i in all_RelationStoichiometry:
			thisType = self._allProductType[self._allProducts[i.reactant_fk_id]]
			self._allRelationStoichiometry[i.id] = { "coeff": float(i.coefficient),
								"location": self._dbLocationId[i.location_fk_id],
								"molecule": self._allProducts[i.reactant_fk_id],
								"form": "mature",
								"type":  thisType
								}

	def _loadPromoters(self):
		self._promoters = []
		self._promoterDbId = {}

		self._checkDatabaseAccess(Promoter)
		all_pr = Promoter.objects.all()
		for i in all_pr:
			self._promoterDbId[i.id] = i.promoter_id
			p = {
				"id":i.promoter_id,
				"name":str(i.name),
				"position":int(i.position),
				"direction":str(i.direction)
			}
			if p["direction"] == "f":
				p["direction"] = '+'
			else:
				p["direction"] = '-'

			self._promoters.append(p)

	def _loadTerminators(self):
		self._terminators = []
		self._terminatorDbId = {}

		self._checkDatabaseAccess(Terminator)
		all_tr = Terminator.objects.all()
		for i in all_tr:
			self._terminatorDbId[i.id] = i.terminator_id
			t = {
				"id":i.terminator_id,
				"name":str(i.name),
				"left":int(i.left),
				"right":int(i.right),
				"rho":str(i.rho_dependent)
			}

			self._terminators.append(t)

	def _loadTranscriptionUnits(self):

		self._transcriptionUnits = []

		#gene
		tu_gene = {}
		self._checkDatabaseAccess(TranscriptionUnitGene)
		all_tg = TranscriptionUnitGene.objects.all()
		for i in all_tg:
			tu = i.transcriptionunit_id_fk_id
			gene = self._geneDbIds[i.gene_id_fk_id]
			if tu in tu_gene:
				tu_gene[tu].append(gene)
			else:
				tu_gene[tu] = [gene]

		#terminator
		tu_tr = {}
		self._checkDatabaseAccess(TranscriptionUnitTerminator)
		all_tt = TranscriptionUnitTerminator.objects.all()
		for i in all_tt:
			tu = i.transcription_unit_id_fk_id
			tr = self._terminatorDbId[i.terminator_id_fk_id]
			if tu in tu_tr:
				tu_tr[tu].append(tr)
			else:
				tu_tr[tu] = [tr]


		self._checkDatabaseAccess(TranscriptionUnit)
		all_tu = TranscriptionUnit.objects.all()
		for i in all_tu:
			t = {
				"id":i.transcription_unit_id,
				"name":str(i.name),
				"left":int(i.left),
				"right":int(i.right),
				"direction":str(i.direction),
				"degradation_rate": float(i.degradation_rate),
				"expression_rate": float(i.expression_rate),
				"promoter_id": self._promoterDbId[i.promoter_id_fk_id],
				"gene_id": tu_gene[i.id],
				"terminator_id": tu_tr[i.id],
			}

			self._transcriptionUnits.append(t)



	def _loadBiomassFractions(self):

		doublingTime = [100, 60, 40, 30, 24]

		self._cellDryMassCompositionData = np.zeros(len(doublingTime),
			dtype = [('doublingTime',				'float64'),
					('proteinMassFraction',			'float64'),
					('rnaMassFraction',				'float64'),
					('dnaMassFraction',				'float64'),
					('lipidMassFraction',			'float64'),
					('lpsMassFraction',				'float64'),
					('mureinMassFraction',			'float64'),
					('glycogenMassFraction',		'float64'),
					('solublePoolMassFraction',		'float64'),
					('inorganicIonMassFraction',	'float64')])

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

		self._cellLipidFractionData = np.zeros(len(lipidIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'float64')])
		self._cellLipidFractionData['metaboliteId'] = lipidIds
		self._cellLipidFractionData['massFraction'] = fracOfLipidMass

		## LPS
		lpsIds = ['colipa[e]']
		fracOfLPSMass = [1.]

		lpsIds = [x[:-3].upper() + x[-3:] for x in lpsIds]
		if abs(sum(fracOfLPSMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellLPSFractionData = np.zeros(len(lpsIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'float64')])
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

		self._cellMureinFractionData = np.zeros(len(mureinIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'float64')])
		self._cellMureinFractionData['metaboliteId'] = mureinIds
		self._cellMureinFractionData['massFraction'] = fracOfMureinMass

		# Glycogen
		glycogenIds = ['glycogen[c]']
		fracOfGlycogenMass = [1.]

		glycogenIds = [x[:-3].upper() + x[-3:] for x in glycogenIds]
		if abs(sum(fracOfGlycogenMass) - 1.0) > 1e-5:
			raise Exception, 'Fractions do not sum to one!\n'

		self._cellGlycogenFractionData = np.zeros(len(glycogenIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'float64')])
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

		self._cellSolublePoolFractionData = np.zeros(len(solublePoolIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'float64')])
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

		self._cellInorganicIonFractionData = np.zeros(len(inorganicIonIds),
			dtype = [('metaboliteId', 'a50'), ('massFraction', 'float64')])
		self._cellInorganicIonFractionData['metaboliteId'] = inorganicIonIds
		self._cellInorganicIonFractionData['massFraction'] = fracInorganicIonMass

	def _loadGenome(self):
		self._translationTable = 11 # E. coli is 11

		all_seq = Chromosome.objects.all()
		genome = ''
		for i in all_seq:
			genome = i.sequence
			break
		self._genomeSeq = genome
		self._genomeLength = len(self._genomeSeq)


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
				g["seq"] = self._genomeSeq[(g["coordinate"]): (g["coordinate"] + g["length"])]
			else:
				g["direction"] = '-'
				g["seq"] = Bio.Seq.Seq(self._genomeSeq[(g["coordinate"] - g["length"] + 1): (g["coordinate"] + 1)]).reverse_complement().tostring()

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

		#uncomment after fixing modified rxn from ecocyc
		'''
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
		'''

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
				"mw" : np.zeros(len(MOLECULAR_WEIGHT_ORDER)),
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
				"mw": np.zeros(len(MOLECULAR_WEIGHT_ORDER))
				}

			if int(i.is_modified):
				if i.id not in rnamodified:
					raise Exception, "%s RNA has no modified form" % i.frame_id_id
				r["modifiedForms"] = rnamodified[i.id]

			if r["type"] == "mRNA":
				r["monomerId"] = self._genes[geneLookup[r["geneId"]]]["monomerId"]

			gene_seq = self._genes[geneLookup[r["geneId"]]]["seq"]
			r["seq"] = Bio.Seq.Seq(gene_seq, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()
			r["ntCount"] = np.array([r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")])
			weight = (
				self._ntWeights["A"] * r["ntCount"][0]
				+ self._ntWeights["C"] * r["ntCount"][1]
				+ self._ntWeights["G"] * r["ntCount"][2]
				+ self._ntWeights["U"] * r["ntCount"][3]
				) + self._rnaEndWeight
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
					"mw": np.zeros(len(MOLECULAR_WEIGHT_ORDER))
				}

				r["seq"] = Bio.Seq.Seq(g["seq"], Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()
				r["ntCount"] = np.array([r["seq"].count("A"), r["seq"].count("C"), r["seq"].count("G"), r["seq"].count("U")])
				weight = (
					self._ntWeights["A"] * r["ntCount"][0]
					+ self._ntWeights["C"] * r["ntCount"][1]
					+ self._ntWeights["G"] * r["ntCount"][2]
					+ self._ntWeights["U"] * r["ntCount"][3]
					) + self._rnaEndWeight
				index = self._whichRna(r['id'], r['type'])
				r["mw"][index] = weight

				self._rnas.append(r)


	def _whichRna(self, rnaId, rnaType):
		if rnaType == 'miscRNA':
			return MOLECULAR_WEIGHT_ORDER['miscRNA']
		if rnaType == 'tRNA':
			return MOLECULAR_WEIGHT_ORDER['tRNA']
		if rnaType == 'mRNA':
			return MOLECULAR_WEIGHT_ORDER['mRNA']
		if rnaType == "rRNA" and rnaId.startswith("RRL"):
			return MOLECULAR_WEIGHT_ORDER['23srRNA']
		if rnaType == "rRNA" and rnaId.startswith("RRS"):
			return MOLECULAR_WEIGHT_ORDER['16srRNA']
		if rnaType == "rRNA" and rnaId.startswith("RRF"):
			return MOLECULAR_WEIGHT_ORDER['5srRNA']

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
				"mw" : np.zeros(len(MOLECULAR_WEIGHT_ORDER)),
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
		#monomermod = self._loadModifiedProteinMonomers()
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
				#"modifiedForms": [],
				"comments": self._allComments[i.comment_fk_id],
				"seq": "",
				"codingRnaSeq" : "",
				"aaCount": np.zeros(21),
				"mw": -1,
				"rnaId": self._genes[geneLookup[gene_frame_id]]["rnaId"]
			}
			#uncomment if need to load modified form
			'''
			if int(i.is_modified):	#TODO: Check after update monomer_modified by Nick
				if i.id in monomermod:
					p["modifiedForms"] = monomermod[i.id]
				else:
					raise Exception, "modified Monomer Absent %s" % p["id"]
					#print p["id"], i.id, i.name

			'''

			if gene_frame_id in genesplices:
				baseSequence = Bio.Seq.Seq("", Bio.Alphabet.IUPAC.IUPACUnambiguousDNA())
				baseSequence = baseSequence + self._genomeSeq[genesplices[gene_frame_id]['start1']-1:genesplices[gene_frame_id]['stop1']] + self._genomeSeq[genesplices[gene_frame_id]['start2']-1:genesplices[gene_frame_id]['stop2']]

				if self._genes[geneLookup[gene_frame_id]]["direction"] == "-":
					baseSequence = baseSequence.reverse_complement()
				baseSequence = baseSequence.tostring()

			else:
				baseSequence = self._genes[geneLookup[gene_frame_id]]['seq']

			p["seq"] = Bio.Seq.Seq(baseSequence, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).translate(table = self._translationTable).tostring()
			p["codingRnaSeq"] = Bio.Seq.Seq(baseSequence, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA()).transcribe().tostring()

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
			p["aaCount"] = np.array([tmp["A"], tmp["R"], tmp["N"], tmp["D"], tmp["C"],
							tmp["E"], tmp["Q"], tmp["G"], tmp["H"], tmp["I"],
							tmp["L"], tmp["K"], tmp["M"], tmp["F"], tmp["P"],
							tmp["U"], tmp["S"], tmp["T"], tmp["W"], tmp["Y"], tmp["V"]
							])

			p["mw"] = self._polypeptideEndWeight
			for aa in p["seq"]: p["mw"] += self._aaWeights[aa]

			self._proteins.append(p)

		#uncomment if need to load modified form
		'''
		##update FK of modified proteins
		for i in self._modifiedProteins:
			i["unmodifiedForm"] = proteinDbIds[i["unmodifiedForm"]]
		'''

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
				"mw" : np.zeros(len(MOLECULAR_WEIGHT_ORDER)),
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
		#complexMod = self._loadModifiedProteinComplexes()
		deletedComplexes = COMPLEXES_REQUIRE_MODIFIED + COMPLEXES_NOT_FORMED
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
			if self._allProducts[i.protein_complex_id] in deletedComplexes: continue #TODO: No need after fixing DB

			complexDbIds[i.id] = self._allProducts[i.protein_complex_id]
			p = {
				"id": self._allProducts[i.protein_complex_id],
				"reactionId": self._allProducts[i.protein_complex_id]+'_RXN',
				"name": i.name,
				#"modifiedForms": [],
				"location": self._dbLocationId[i.location_fk_id],
				"mw": np.zeros(len(MOLECULAR_WEIGHT_ORDER)),
				"comments": self._allComments[i.comment_fk_id]
			}
			#uncomment if need to load modified form
			'''
			if i.modified_form:
				if i.id in complexMod:
					p["modifiedForms"] = complexMod[i.id]
				else:
					raise Exception, "modification form absent for complex"
			'''
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

		#uncomment if need to load modified form
		'''
		##update FK of modified complexes
		for i in self._modifiedComplexes:
			i["unmodifiedForm"] = complexDbIds[i["unmodifiedForm"]]
		'''

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
					"lb": float(i.lower_bound),
					"kcat":i.kcat if i.kcat != -1 else None
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

	'''
	def _loadTranscriptionUnits(self):

		self._transcriptionUnits = []

		#gene
		tu_gene = {}
		self._checkDatabaseAccess(TranscriptionUnitGene)
		all_tg = TranscriptionUnitGene.objects.all()
		for i in all_tg:
			tu = i.transcriptionunit_id_fk_id
			gene = self._geneDbIds[i.gene_id_fk_id]
			if tu in tu_gene:
				tu_gene[tu].append(gene)
			else:
				tu_gene[tu] = [gene]


		tu_pr = {}

		self._checkDatabaseAccess(TranscriptionUnit)
		all_tu = TranscriptionUnit.objects.all()
		for i in all_tu:
			t = {
				"id":i.transcription_unit_id,
				"name":str(i.name),
				"left":int(i.left),
				"right":int(i.right),
				"direction":str(i.direction),
				"degradation_rate": float(i.degradation_rate),
				"expression_rate": float(i.expression_rate),
				"promoter_id": self._promoterDbId[i.promoter_id_fk_id],
				"gene_id": tu_gene[i.id]
			}
			if t["promoter_id"] in tu_pr:
				tu_pr[t["promoter_id"]].append(t["id"])
			else:
				tu_pr[t["promoter_id"]] = [t["id"]]

			self._transcriptionUnits.append(t)

		#Add TU info in promoters
		for p in self._promoters:
			p['TU'] = tu_pr[p['id']]

	def _countATinPromoters(self):

		geneLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._genes)])
		tuLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._transcriptionUnits)])

		#calculate AT counts for each genes associated with each promoters
		genes_pr = {}
		for g in self._genes:
			genes_pr[g['id']] = []

		for p in self._promoters:
			p['TA_count'] = (p['seq'].count('A') + p['seq'].count('T'))/float(len(p['seq'])) *100
			allgenes = []
			for t in p['TU']:
				genes = self._transcriptionUnits[tuLookUp[t]]['gene_id']
				for g in genes:
					if g in allgenes:
						continue
					allgenes.append(g)
					frame_id = self._genes[geneLookUp[g]]['id']
					genes_pr[frame_id].append(p['TA_count'])
		total = 0
		for g in self._genes:
			if len(genes_pr[g['id']]) == 0:
				#print g['id']
				total = total + 1
				continue
			x = sum(genes_pr[g['id']])/float(len(genes_pr[g['id']]))
			#print g['id'],'\t',g['name'],'\t',g['symbol'],'\t',x ,'\t', len(genes_pr[g['id']])
		#print total
	'''

	def _loadMetaboliteConcentrations(self):
		# TODO: move data to SQL and load here

		self._metaboliteConcentrations = [
			(metaboliteID.upper(), concentration)
			for metaboliteID, concentration in METABOLITE_CONCENTRATIONS.viewitems()
			]


	def _calcMolecularWeightFromRxn(self):

		complexReactionLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._complexationReactions)])
		modificationReactionLookUp = dict([(x[1]["id"], x[0]) for x in enumerate(self._modificationReactions)])
		met = dict([(x["id"], x['mw7.2']) for x in self._metabolites])
		rna = dict([(x["id"], x['mw']) for x in self._rnas])
		monomer = dict([(x["id"], x['mw']) for x in self._proteins])

		products = {}

		for i in self._proteinComplexes:
			products[i['id']] = self._complexationReactions[complexReactionLookUp[i['reactionId']]]['stoichiometry']
		for i in self._modifiedRnas:
			r = i['reactionId'][len(i['reactionId'])-1] #last reaction; considering only 1 reaction
			products[i['id']] = self._modificationReactions[modificationReactionLookUp[r]]['stoichiometry']

		#initialize matrixes
		counts = {} # number of uncomputable products in the rxn
		dependant = {} # [i][j] where product j depends on product i
		queue = [] # list of products which have no dependancy

		for p in products:
			dependant[p] = []

		for p in products:
			tmpCount = 0
			for r in products[p]: #for a molecule in the reaction
				if r['molecule'] == p: continue
				if r['molecule'] in products:
					tmpCount = tmpCount + 1
					if p not in dependant[r['molecule']]:
						dependant[r['molecule']].append(p)

				elif r['type'] != 'metabolite' and r['type'] != 'rna' and r['type'] != 'proteinmonomers':
					#print '\''+frame_id+'\',',i['molecule'], i['type']
					raise Exception, "%s unknown molecule while calculating MW" % r['molecule']

			counts[p] = tmpCount
			if tmpCount == 0:
				queue.append(p) # can be computable

		#topological order
		totalProducts = 0
		mw = {}
		while (len(queue)):
			newQueue = []
			for i in queue:
				#calculate weight
				weight = np.zeros(len(MOLECULAR_WEIGHT_ORDER))
				sign_wt = 1
				for m in products[i]:
					if m['molecule'] in products:
						if i == m['molecule']:
							sign_wt = m['coeff']

						else:
							weight += mw[m['molecule']] * m['coeff']

					elif m['type'] == 'metabolite':
						index = MOLECULAR_WEIGHT_ORDER['metabolite']
						weight[index] += met[m['molecule'].upper()] * m['coeff']

					elif m['type'] == 'rna':
						weight += rna[m['molecule']] * m['coeff']

					elif m['type'] == 'proteinmonomers':
						index = MOLECULAR_WEIGHT_ORDER['protein']
						weight[index] += monomer[m['molecule']] * m['coeff']

					else:
					 	raise Exception, "%s dependant molecule while calculating MW" % m['molecule']

				mw[i] = weight/ (sign_wt * (-1))

				#update dependant and count
				for j in dependant[i]:
					counts[j] = counts[j] -1
					if counts[j] == 0:
						newQueue.append(j)
			#update queue
			totalProducts = totalProducts + len(queue)
			queue = newQueue

 		if totalProducts < len(products):
			print 'there are cycles in rxn!', totalProducts,len(products),len(mw)

		#update the MW for proteinComplexes and modifiedRnas
		for i in self._modifiedRnas:
			i['mw'] = mw[i['id']]

		for i in self._proteinComplexes:
			i['mw'] = mw[i['id']]


	def _calcKCat(self, enzId, vMax, units):
		if enzId == None or vMax == None:
			return np.NaN

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
		self._constantData = {}
		self._constantData['nAvogadro'] = scipy.constants.Avogadro * 1 / units.mol

	def _loadParameters(self):
		self._parameterData = {}
		self._parameterData['cellCycleLen'] = 3600*units.s
		self._parameterData['avgCellDryMass'] = 258*units.fg
		self._parameterData['rnaPolymeraseElongationRate'] = 42*units.nt/units.s
		self._parameterData['rnaPolymeraseElongationRateFast'] = 80*units.nt/units.s
		self._parameterData['ribosomeElongationRate'] = 16*units.aa/units.s
		self._parameterData['fracInitFreeNTPs'] = 0.0015
		self._parameterData['fracInitFreeAAs'] = 0.001
		self._parameterData['avgCellCellCycleProgress'] = 0.44
		self._parameterData['timeStep'] = 1*units.s

	def _loadComputeParameters(self):
		self._parameterData['avgCellToInitalCellConvFactor'] = np.exp(np.log(2) * self._parameterData['avgCellCellCycleProgress'])
		self._parameterData['avgCellDryMassInit'] = self._parameterData['avgCellDryMass'] / self._parameterData['avgCellToInitalCellConvFactor']
		self._parameterData['avgCellWaterMass'] = (self._parameterData['avgCellDryMass'] / self._parameterData['cellDryMassFraction']) * self._parameterData['cellWaterMassFraction']
		self._parameterData['avgCellWaterMassInit'] = self._parameterData['avgCellWaterMass'] / self._parameterData['avgCellToInitalCellConvFactor']


	## -- Build functions -- ##

	def _buildSequence(self):
		self.genomeSeq = self._genomeSeq
		self.genomeLength = self._genomeLength
		self.genome_A_count = self.genomeSeq.count("A")
		self.genome_T_count = self.genomeSeq.count("T")
		self.genome_G_count = self.genomeSeq.count("G")
		self.genome_C_count = self.genomeSeq.count("C")


	def _buildCompartments(self):
		self._compartmentData = np.zeros(len(self._compartmentList),
			dtype = [('compartmentId','a20'),('compartmentAbbreviation', 'a1')])

		# Load data into structured array
		self._compartmentData['compartmentId']				= [x['id'] for x in self._compartmentList]
		self._compartmentData['compartmentAbbreviation']	= [x['abbrev'] for x in self._compartmentList]
		self.compartments = self._compartmentData
		self.nCompartments 	= len(self._compartmentList)


	def _buildBulkMolecules(self):
		# TODO: modularize this logic

		size = (
			len(self._metabolites)*len(self._compartmentList)
			+ len(self._rnas)
			+ len(self._proteins)
			+ len(self._proteinComplexes)
			+ len(self._polymerized)*len(self._compartmentList)
			)

		bulkMolecules = np.zeros(
			size,
			dtype = [
				("moleculeId", "a50"),
				('compartment',	 "a1"),
				("mass", "{}f8".format(len(MOLECULAR_WEIGHT_ORDER))),
				]
			)

		# Set metabolites
		lastMetaboliteIdx = len(self._metabolites) * len(self._compartmentList)

		compartmentAbbreviations = [compartment['abbrev'] for compartment in self._compartmentList]
		metaboliteIds = [metabolite['id'] for metabolite in self._metabolites]

		bulkMolecules['moleculeId'][0:lastMetaboliteIdx] = [
			'{}[{}]'.format(metaboliteId, compartmentAbbreviation)
			for compartmentAbbreviation in compartmentAbbreviations
			for metaboliteId in metaboliteIds
			]

		metaboliteMassIdxs = np.empty(lastMetaboliteIdx, np.int64)

		metaboliteMassIdxs.fill(MOLECULAR_WEIGHT_ORDER["metabolite"])

		for index, metabolite in enumerate(bulkMolecules[0:lastMetaboliteIdx]):
			if metabolite["moleculeId"].startswith("H2O["):
				metaboliteMassIdxs[index] = MOLECULAR_WEIGHT_ORDER["water"]

		bulkMolecules['mass'][np.arange(lastMetaboliteIdx), metaboliteMassIdxs] = [
			metabolite['mw7.2']
			for compartmentIndex in range(len(self._compartmentList))
			for metabolite in self._metabolites
			]

		# Set RNA
		lastRnaIdx = len(self._rnas) + lastMetaboliteIdx

		bulkMolecules['moleculeId'][lastMetaboliteIdx:lastRnaIdx] = [
			'{}[{}]'.format(rna['id'], rna['location']) for rna in self._rnas
			]

		bulkMolecules['mass'][lastMetaboliteIdx:lastRnaIdx, :] = [
			rna['mw'] for rna in self._rnas
			]

		# Set proteins
		lastProteinMonomerIdx = len(self._proteins) + lastRnaIdx

		bulkMolecules['moleculeId'][lastRnaIdx:lastProteinMonomerIdx] = [
			'{}[{}]'.format(protein['id'], protein['location'])
			for protein in self._proteins
			]

		bulkMolecules['mass'][lastRnaIdx:lastProteinMonomerIdx, MOLECULAR_WEIGHT_ORDER["protein"]] = [
			protein['mw'] for protein in self._proteins
			]

		# Set complexes
		lastComplexIdx = len(self._proteinComplexes) + lastProteinMonomerIdx

		bulkMolecules['moleculeId'][lastProteinMonomerIdx:lastComplexIdx] = [
			'{}[{}]'.format(complex_['id'],complex_['location']) for complex_ in self._proteinComplexes
			]

		bulkMolecules['mass'][lastProteinMonomerIdx:lastComplexIdx, :] = [
			complex_['mw'] for complex_ in self._proteinComplexes
			]

		# Set polymerized

		lastPolymerizedIndex = len(self._polymerized)*len(self._compartmentList) + lastComplexIdx

		polymerizedIDs = [entry["id"] for entry in self._polymerized]

		bulkMolecules["moleculeId"][lastComplexIdx:lastPolymerizedIndex] = [
			'{}[{}]'.format(polymerizedID, compartmentAbbreviation)
			for compartmentAbbreviation in compartmentAbbreviations
			for polymerizedID in polymerizedIDs
			]

		masses = [
			entry["mw"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in self._polymerized
			]

		massIndexes = [
			entry["mass key"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in self._polymerized
			]

		bulkMolecules["mass"][range(lastComplexIdx, lastPolymerizedIndex), massIndexes] = masses
		# NOTE: the use of range above is intentional

		# Add units to values
		field_units = {
			"moleculeId"		:	None,
			"mass"				:	units.g / units.mol,
			'compartment'		:	None,
			}

		self.bulkMolecules = UnitStructArray(bulkMolecules, field_units)


	def _buildGeneData(self):
		self.geneData = np.zeros(len(self._genes),
			dtype = [('name'				,	'a50'),
					#('coordinate'			,	'int64'),
					#('length'				,	'int64'),
					#('positiveDirection'	,	'bool'),
					('rnaId'                ,   'a50'),
					('endCoordinate'		,	'int64')])

		self.geneData['name'] = [x['id'] for x in self._genes]
		self.geneData['rnaId'] = [x['rnaId'] for x in self._genes]
		#self.geneData['coordinate'] = [x['coordinate'] for x in self._genes]
		#self.geneData['length'] = [x['length'] for x in self._genes]
		#self.geneData['positiveDirection'] = [True if x['direction'] == '+' else False for x in self._genes]
		self.geneData['endCoordinate'] = [(x['coordinate'] + x['length']) % self.genomeLength if x['direction'] == '+' else (x['coordinate'] - x['length']) % self.genomeLength for x in self._genes]


	def _buildBulkChromosome(self):
		count_dnaAbox_at_oriC = 5
		count_dnaABindingSites_at_oriC = 5

		dnaA_mass = [x['mw'] for x in self._proteins if x['id'] == 'PD03831'][0]
		atp_mass = [x['mw7.2'] for x in self._metabolites if x['id'] == 'ATP'][0]
		adp_mass = [x['mw7.2'] for x in self._metabolites if x['id'] == 'ADP'][0]
		dnaA_atp_mass = dnaA_mass + atp_mass
		dnaA_adp_mass = dnaA_mass + adp_mass

		size = len(self._genes) + count_dnaABindingSites_at_oriC + count_dnaABindingSites_at_oriC + count_dnaABindingSites_at_oriC
		bulkChromosome = np.zeros(size,
			dtype = [("moleculeId", 			"a50"),
					('compartment',				"a1"),
					("mass", "{}f8".format(len(MOLECULAR_WEIGHT_ORDER))),
					("isGene",					"bool"),
					("isDnaABox",				"bool"),
					("isDnaABox_atp_polymer",	"bool"),
					("isDnaABox_adp_polymer",	"bool")]
					)

		# Set genes
		lastGeneIdx = len(self._genes)
		bulkChromosome['moleculeId'][0:lastGeneIdx] = [x['id'] for x in self._genes]
		bulkChromosome['compartment'][0:lastGeneIdx] = 'n'
		bulkChromosome['isGene'][0:lastGeneIdx] = True

		# Set dnaA box
		lastDnaAIdx = lastGeneIdx + count_dnaAbox_at_oriC
		bulkChromosome['compartment'][lastGeneIdx:lastDnaAIdx] = 'n'
		bulkChromosome['moleculeId'][lastGeneIdx:lastDnaAIdx] = ['R1_dnaA',
																'R2_dnaA',
																'R3_dnaA',
																'R4_dnaA',
																'R5_dnaA']
		bulkChromosome['isDnaABox'][lastGeneIdx:lastDnaAIdx] = True

		# Set dnaA box dnaA-ATP polymer
		lastDnaAATPSiteIdx = lastDnaAIdx + count_dnaABindingSites_at_oriC
		bulkChromosome['compartment'][lastDnaAIdx:lastDnaAATPSiteIdx] = 'n'
		bulkChromosome['moleculeId'][lastDnaAIdx:lastDnaAATPSiteIdx] = ['R1_dnaA_atp_polymer',
																'R2_dnaA_atp_polymer',
																'R3_dnaA_atp_polymer',
																'R4_dnaA_atp_polymer',
																'R5_dnaA_atp_polymer']
		bulkChromosome['mass'][lastDnaAIdx:lastDnaAATPSiteIdx, MOLECULAR_WEIGHT_ORDER["protein"]] = dnaA_atp_mass
		bulkChromosome['isDnaABox_atp_polymer'][lastDnaAIdx:lastDnaAATPSiteIdx] = True

		# Set dnaA box dnaA-ADP polymer
		lastDnaAADPSiteIdx = lastDnaAATPSiteIdx + count_dnaABindingSites_at_oriC
		bulkChromosome['compartment'][lastDnaAATPSiteIdx:lastDnaAADPSiteIdx] = 'n'
		bulkChromosome['moleculeId'][lastDnaAATPSiteIdx:lastDnaAADPSiteIdx] = ['R1_dnaA_adp_polymer',
																'R2_dnaA_adp_polymer',
																'R3_dnaA_adp_polymer',
																'R4_dnaA_adp_polymer',
																'R5_dnaA_adp_polymer']
		bulkChromosome['mass'][lastDnaAATPSiteIdx:lastDnaAADPSiteIdx, MOLECULAR_WEIGHT_ORDER["protein"]] = dnaA_adp_mass
		bulkChromosome['isDnaABox_adp_polymer'][lastDnaAATPSiteIdx:lastDnaAADPSiteIdx] = True


		# Add units to values
		field_units = {
			"moleculeId"			:	None,
			"mass"					:	units.g / units.mol,
			'compartment'			:	None,
			'isGene'				:	None,
			"isDnaABox"				:	None,
			"isDnaABox_atp_polymer"	:	None,
			"isDnaABox_adp_polymer"	:	None,
			}
		self.bulkChromosome = UnitStructArray(bulkChromosome, field_units)

	def _buildUniqueMolecules(self):

		self.uniqueMoleculeDefinitions = collections.OrderedDict([
			("activeRnaPoly", {
				'rnaIndex' : 'i8',
				'transcriptLength' : 'i8'
				}),
			("activeRibosome", {
				'proteinIndex' : 'i8',
				'peptideLength': 'i8'
				}),

			("dnaPolymerase", {
				'chromosomeLocation' : 'i8',
				'directionIsPositive' : 'bool',
				'isLeading' : 'bool'
				}),
			])

		rnaPolyComplexMass = self.bulkMolecules["mass"][self.bulkMolecules["moleculeId"] == "APORNAP-CPLX[c]"].asNumber()

		# TODO: This is a bad hack that works because in the fitter
		# I have forced expression to be these subunits only
		ribosomeSubunits = [self.s30_fullComplex, self.s50_fullComplex]

		ribosomeMass = sum(
			entry["mass"] for entry in self.bulkMolecules.struct_array
			if entry["moleculeId"] in ribosomeSubunits
			)

		dnaPolyMass = np.zeros_like(rnaPolyComplexMass) # NOTE: dnaPolymerases currently have no mass

		uniqueMoleculeMasses = np.zeros(
			shape = len(self.uniqueMoleculeDefinitions),
			dtype = [
				('moleculeId', 'a50'),
				("mass", "{}f8".format(len(MOLECULAR_WEIGHT_ORDER))),
				]
			)

		uniqueMoleculeMasses["moleculeId"] = self.uniqueMoleculeDefinitions.keys()

		uniqueMoleculeMasses["mass"] = np.vstack([
			rnaPolyComplexMass,
			ribosomeMass,
			dnaPolyMass
			])

		self.uniqueMoleculeMasses = UnitStructArray(
			uniqueMoleculeMasses,
			{"moleculeId":None, "mass":units.g / units.mol}
			)

		# TODO: add the ability to "register" a bulk molecule as a unique
		# molecule to handle most of the above logic


	def _buildRnaExpression(self):
		normalizedRnaExpression = np.zeros(sum(1 for x in self._rnas),
			dtype = [('rnaId',		'a50'),
					('expression',	'float64'),
					('isMRna',		'bool'),
					('isMiscRna',	'bool'),
					('isRRna',		'bool'),
					('isTRna',		'bool'),
					('isRRna23S',	'bool'),
					('isRRna16S',	'bool'),
					('isRRna5S',	'bool')])

		normalizedRnaExpression['rnaId'] 		= ['{}[{}]'.format(x['id'], x['location']) for x in self._rnas]
		normalizedRnaExpression['expression']	= [x['expression'] for x in self._rnas]
		normalizedRnaExpression['expression']	= normalizedRnaExpression['expression'] / np.sum(normalizedRnaExpression['expression'])
		normalizedRnaExpression['isMRna'] = [rna["type"] == "mRNA" for rna in self._rnas]
		normalizedRnaExpression['isMiscRna'] = [rna["type"] == "miscRNA" for rna in self._rnas]
		normalizedRnaExpression['isRRna'] = [rna["type"] == "rRNA" for rna in self._rnas]
		normalizedRnaExpression['isTRna'] = [rna["type"] == "tRNA" for rna in self._rnas]

		self.rnaExpression = UnitStructArray(normalizedRnaExpression,
			{
			'rnaId'		:	None,
			'expression':	None,
			'isMRna'	:	None,
			'isMiscRna'	:	None,
			'isRRna'	:	None,
			'isTRna'	:	None,
			'isRRna23S'	:	None,
			'isRRna16S'	:	None,
			'isRRna5S'	:	None
			})


	def _buildBiomass(self):
		self._coreBiomassData = np.zeros(sum(len(x['biomassInfo']['core']) for x in self._metabolites if len(x['biomassInfo']['core'])),
			dtype = [('metaboliteId', 'a50'),
					('biomassFlux', 	'float64')])

		self._wildtypeBiomassData = np.zeros(sum(len(x['biomassInfo']['wildtype']) for x in self._metabolites if len(x['biomassInfo']['wildtype'])),
			dtype = [('metaboliteId', 'a50'),
					('biomassFlux',		'float64')])

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

		field_units = {'metaboliteId' : None,
				'biomassFlux' : units.mmol / units.g}
		self.coreBiomass 		= UnitStructArray(self._coreBiomassData, field_units)
		self.wildtypeBiomass 	= UnitStructArray(self._wildtypeBiomassData, field_units)

	def _buildBiomassFractions(self):
		field_units = {
			'doublingTime' : units.min,
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

		self.cellDryMassComposition = UnitStructArray(self._cellDryMassCompositionData, field_units)
		self.cellLipidFractionData = self._cellLipidFractionData
		self.cellLPSFractionData = self._cellLPSFractionData
		self.cellMureinFractionData = self._cellMureinFractionData
		self.cellGlycogenFractionData = self._cellGlycogenFractionData
		self.cellSolublePoolFractionData = self._cellSolublePoolFractionData
		self.cellInorganicIonFractionData = self._cellInorganicIonFractionData

	def _buildRnaData(self):
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location']) for rna in self._rnas]

		rnaDegRates = np.log(2) / np.array([rna['halfLife'] for rna in self._rnas]) # TODO: units

		rnaLens = np.array([len(rna['seq']) for rna in self._rnas])

		ntCounts = np.array([
			(rna['seq'].count('A'), rna['seq'].count('C'),
				rna['seq'].count('G'), rna['seq'].count('U'))
			for rna in self._rnas
			])

		expression = np.array([rna['expression'] for rna in self._rnas])

		synthProb = expression * (
			np.log(2) / self._parameterData['cellCycleLen'].asNumber(units.s)
			+ rnaDegRates
			)

		synthProb /= synthProb.sum()

		mws = np.array([rna['mw'] for rna in self._rnas])

		geneIds = np.array([rna['geneId'] for rna in self._rnas])

		size = len(rnaIds)

		is23S = np.zeros(size, dtype = np.bool)
		is16S = np.zeros(size, dtype = np.bool)
		is5S = np.zeros(size, dtype = np.bool)

		for rnaIndex, rna in enumerate(self._rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is23S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is16S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is5S[rnaIndex] = True

		sequences = [rna['seq'] for rna in self._rnas]

		maxSequenceLength = max(len(sequence) for sequence in sequences)

		# TODO: Add units
		self.rnaData = np.zeros(
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
				('isRRna5S', 'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('geneId', 'a50')
				]
			)

		self.rnaData['id'] = rnaIds
		self.rnaData['synthProb'] = synthProb
		self.rnaData['degRate'] = rnaDegRates
		self.rnaData['length'] = rnaLens
		self.rnaData['countsACGU'] = ntCounts
		self.rnaData['mw'] = mws.sum(axis = 1)
		self.rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in self._rnas]
		self.rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in self._rnas]
		self.rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in self._rnas]
		self.rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in self._rnas]
		self.rnaData['isRRna23S'] = is23S
		self.rnaData['isRRna16S'] = is16S
		self.rnaData['isRRna5S'] = is5S
		self.rnaData['sequence'] = sequences
		self.rnaData['geneId'] = geneIds

		field_units = {
			'id'		:	None,
			'synthProb' :	None,
			'degRate'	:	1 / units.s,
			'length'	:	units.nt,
			'countsACGU':	units.nt,
			'mw'		:	units.g / units.mol,
			'isMRna'	:	None,
			'isMiscRna'	:	None,
			'isRRna'	:	None,
			'isTRna'	:	None,
			'isRRna23S'	:	None,
			'isRRna16S'	:	None,
			'isRRna5S'	:	None,
			'sequence'  :   None,
			'geneId'	:	None,
			}


		self.rnaData = UnitStructArray(self.rnaData, field_units)
		self.getTrnaAbundanceData = getTrnaAbundanceAtGrowthRate

	def _buildMonomerData(self):
		ids = ['{}[{}]'.format(protein['id'], protein['location'])
			for protein in self._proteins]

		rnaIds = []

		for protein in self._proteins:
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
		sequences = []

		for protein in self._proteins:
			sequence = protein['seq']

			counts = []

			for aa in self._aaWeights.viewkeys(): # TODO: better way to get AA ids?
				counts.append(
					sequence.count(aa)
					)

			lengths.append(len(sequence))
			aaCounts.append(counts)
			sequences.append(sequence)

		maxSequenceLength = max(len(seq) for seq in sequences)

		mws = np.array([protein['mw'] for protein in self._proteins])

		size = len(rnaIds)

		nAAs = len(aaCounts[0])

		# Calculate degradation rates based on N-rule
		# TODO: citation
		fastRate = (np.log(2) / (2*units.min)).asUnit(1 / units.s)
		slowRate = (np.log(2) / (10*60*units.min)).asUnit(1 / units.s)

		fastAAs = ["R", "K", "F", "L", "W", "Y"]
		slowAAs = ["H", "I", "D", "E", "N", "Q", "C", "A", "S", "T", "G", "V", "M"]
		noDataAAs = ["P", "U"]

		NruleDegRate = {}
		NruleDegRate.update(
			(fastAA, fastRate) for fastAA in fastAAs
			)
		NruleDegRate.update(
			(slowAA, slowRate) for slowAA in slowAAs
			)
		NruleDegRate.update(
			(noDataAA, slowRate) for noDataAA in noDataAAs
			) # Assumed slow rate because of no data

		# Build list of ribosomal proteins
		# Give all ribosomal proteins the slowAA rule
		ribosomalProteins = []
		ribosomalProteins.extend([x[:-3] for x in S30_ALL_PROTEINS])
		ribosomalProteins.extend([x[:-3] for x in S50_ALL_PROTEINS])

		degRate = np.zeros(len(self._proteins))
		for i,m in enumerate(self._proteins):
			if m['id'] not in ribosomalProteins:
				degRate[i] = NruleDegRate[m['seq'][0]].asNumber()
			else:
				degRate[i] = slowRate.asNumber()

		self.monomerData = np.zeros(
			size,
			dtype = [
				('id', 'a50'),
				('rnaId', 'a50'),
				('degRate', 'f8'),
				('length', 'i8'),
				('aaCounts', '{}i8'.format(nAAs)),
				('mw', 'f8'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				]
			)

		self.monomerData['id'] = ids
		self.monomerData['rnaId'] = rnaIds
		self.monomerData['degRate'] = degRate
		self.monomerData['length'] = lengths
		self.monomerData['aaCounts'] = aaCounts
		self.monomerData['mw'] = mws
		self.monomerData['sequence'] = sequences

		self.aaIDs = AMINO_ACID_1_TO_3_ORDERED.values()
		self.aaIDs_singleLetter = AMINO_ACID_1_TO_3_ORDERED.keys()

		field_units = {
			'id'		:	None,
			'rnaId'		:	None,
			'degRate'	:	1 / units.s,
			'length'	:	units.aa,
			'aaCounts'	:	units.aa,
			'mw'		:	units.g / units.mol,
			'sequence'  :   None
			}

		self.monomerData = UnitStructArray(self.monomerData, field_units)


	def _buildRnaIndexToMonomerMapping(self):
		self.rnaIndexToMonomerMapping = np.array([np.where(x == self.rnaData["id"])[0][0] for x in self.monomerData["rnaId"]])


	def _buildMonomerIndexToRnaMapping(self):
		self.monomerIndexToRnaMapping = np.array([np.where(x == self.monomerData["rnaId"])[0][0] for x in self.rnaData["id"] if len(np.where(x == self.monomerData["rnaId"])[0])])


	def _buildRnaIndexToGeneMapping(self):
		self.rnaIndexToGeneMapping = np.array([np.where(x + "[c]" == self.rnaData["id"])[0][0] for x in self.geneData["rnaId"]])


	def _buildRibosomeData(self):
		self.s30_proteins = S30_PROTEINS
		self.s30_16sRRNA = S30_16S_RRNAS
		self.s30_fullComplex = S30_FULLCOMPLEX
		self.s50_proteins = S50_PROTEINS
		self.s50_proteinComplexes = S50_PROTEIN_COMPLEXES
		self.s50_23sRRNA = S50_23S_RRNAS
		self.s50_5sRRNA = S50_5S_RRNAS
		self.s50_fullComplex = S50_FULLCOMPLEX


	def _buildMetabolitePools(self):
		CELL_DENSITY = 1.1e3 # g/L
		# from Baldwin WW, Myer R, Powell N, Anderson E, Koch AL. Buoyant
		# density of Escherichia coli is determined solely by the osmolarity
		# of the culture medium. Arch Microbiol. 1995 Aug164(2):155-7 p.156
		# fig.3 & fig.2 (retrieved from Bionumbers)

		# Create vector of metabolite pools (concentrations)

		# Since the data only covers certain metabolites, we need to rationally
		# expand the dataset to include the other molecules in the biomass
		# function.

		# First, load in metabolites that do have concentrations, then assign
		# compartments according to those given in the biomass objective.  Or,
		# if there is no compartment, assign it to the cytoplasm.

		metaboliteIDs = []
		metaboliteConcentrations = []

		wildtypeIDs = self._wildtypeBiomassData["metaboliteId"].tolist()

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			}

		for metaboliteID, concentration in self._metaboliteConcentrations:
			if metaboliteID in wildtypeIDtoCompartment:
				metaboliteIDs.append(
					metaboliteID + wildtypeIDtoCompartment[metaboliteID]
					)

			elif metaboliteID == "23CAMP":
				metaboliteIDs.append(
					metaboliteID + "[p]"
					)

			else:
				metaboliteIDs.append(
					metaboliteID + "[c]"
					)

			metaboliteConcentrations.append(concentration)

		# Calculate the following assuming 60 min doubling time

		initWaterMass = self.avgCellWaterMassInit.asNumber(units.g)
		initDryMass = self.avgCellDryMassInit.asNumber(units.g)

		initCellMass = initWaterMass + initDryMass

		initCellVolume = initCellMass / CELL_DENSITY # L

		massFractions = self.cellDryMassComposition[
			self.cellDryMassComposition["doublingTime"].asNumber(units.min) == 60.0
			].fullArray()

		for entry in self.cellGlycogenFractionData:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["glycogenMassFraction"][0]
			molWeight = self.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self.cellMureinFractionData:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["mureinMassFraction"][0]
			molWeight = self.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self.cellLPSFractionData:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["lpsMassFraction"][0]
			molWeight = self.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self.cellLipidFractionData:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["lipidMassFraction"][0]
			molWeight = self.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self.cellInorganicIonFractionData:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["inorganicIonMassFraction"][0]
			molWeight = self.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self.cellSolublePoolFractionData:
			metaboliteID = entry["metaboliteId"]

			if metaboliteID not in metaboliteIDs:
				massFrac = entry["massFraction"] * massFractions["solublePoolMassFraction"][0]
				molWeight = self.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

				massInit = massFrac * initDryMass
				molesInit = massInit/molWeight

				concentration = molesInit / initCellVolume

				metaboliteIDs.append(metaboliteID)
				metaboliteConcentrations.append(concentration)


		# ILE/LEU: split reported concentration according to their relative abundances

		aaAbundances = self.monomerData["aaCounts"].asNumber().sum(axis = 0)
		# TODO: more thorough estimate of abundance or some external data point (neidhardt?)

		ileAbundance = aaAbundances[self.aaIDs.index("ILE-L[c]")]
		leuAbundance = aaAbundances[self.aaIDs.index("LEU-L[c]")]

		ILE_LEU_CONCENTRATION = 3.0e-4 # mmol/L

		ileRelative = ileAbundance / (ileAbundance + leuAbundance)
		leuRelative = 1 - ileRelative

		metaboliteIDs.append("ILE-L[c]")
		metaboliteConcentrations.append(ileRelative * ILE_LEU_CONCENTRATION)

		metaboliteIDs.append("LEU-L[c]")
		metaboliteConcentrations.append(leuRelative * ILE_LEU_CONCENTRATION)

		# CYS/SEC/GLY: fit a relative abundance:concentration line (L1 norm)
		# with other amino acids and solve for these

		aaConcentrations = []
		# aaAbundancesWithConcentrations = []

		for aaIndex, aaID in enumerate(self.aaIDs):
			if aaID in metaboliteIDs:
				metIndex = metaboliteIDs.index(aaID)
				aaConcentrations.append(metaboliteConcentrations[metIndex])
				# aaAbundancesWithConcentrations.append(aaAbundances[aaIndex])

		# TODO: implement L1-norm minimization

		# for now: just choosing and assigning the smallest value

		aaSmallestConc = min(aaConcentrations)

		# HACK: min conc. doesn't work here
		metaboliteIDs.append("GLY[c]")
		metaboliteConcentrations.append(
			metaboliteConcentrations[metaboliteIDs.index("ALA-L[c]")]
			)

		metaboliteIDs.append("CYS-L[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		metaboliteIDs.append("SEC-L[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		# DGTP: set to smallest of all other DNTP concentrations

		dntpConcentrations = []
		# dntpAbundancesWithConcentrations = []

		for dntpIndex, dntpID in enumerate(self.dNtpIds):
			if dntpID in metaboliteIDs:
				metIndex = metaboliteIDs.index(dntpID)
				dntpConcentrations.append(metaboliteConcentrations[metIndex])

		dntpSmallestConc = min(dntpConcentrations)

		metaboliteIDs.append("DGTP[c]")
		metaboliteConcentrations.append(dntpSmallestConc)

		# H2O: reported water content of E. coli

		h2oMolWeight = self.getMass(["H2O[c]"])[0].asNumber(units.g / units.mol)
		h2oMoles = initWaterMass/h2oMolWeight

		h2oConcentration = h2oMoles / initCellVolume

		metaboliteIDs.append("H2O[c]")
		metaboliteConcentrations.append(h2oConcentration)

		# H: reported pH

		ECOLI_PH = 7.2

		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append("H[c]")
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI: multiple sources report 0.5 mM

		PPI_CONCENTRATION = 0.5e-3 # M, multiple sources

		# NOTE: Nick says that the physiological levels of PPI are very low - investigate this

		metaboliteIDs.append("PPI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		unaccounted = set(wildtypeIDs) - set(metaboliteIDs)

		assert len(unaccounted) == 0

		# Add byproducts with no annotated concentration to force recycling

		metaboliteIDs.append("UMP[c]")
		metaboliteConcentrations.append(2.40e-5)

		# Other quantities to consider:
		# - (d)NTP byproducts not currently included

		self.metabolitePoolIDs = metaboliteIDs
		self.metabolitePoolConcentrations = units.mol/units.L * np.array(metaboliteConcentrations)
		self.cellDensity = units.g/units.L * CELL_DENSITY


	def _buildTranscription(self):
		from wholecell.utils.polymerize import PAD_VALUE

		sequences = self.rnaData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.rnaData["length"].asNumber().max()
			+ self.rnaPolymeraseElongationRate.asNumber(units.nt / units.s)
			)

		self.transcriptionSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.transcriptionSequences.fill(PAD_VALUE)

		ntMapping = {ntpId:i for i, ntpId in enumerate(["A", "C", "G", "U"])}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.transcriptionSequences[i, j] = ntMapping[letter]

		# TODO: (URGENT) unify peptide weight calculations!
		self.transcriptionMonomerWeights = (
			(
				self.getMass(self.ntpIds)
				- self.getMass(["PPI[c]"])
				)
			/ self.nAvogadro
			).asNumber(units.fg)

		self.transcriptionEndWeight = (self.getMass(["PPI[c]"]) / self.nAvogadro).asNumber(units.fg)


	def _buildTranslation(self):
		from wholecell.utils.polymerize import PAD_VALUE

		sequences = self.monomerData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.monomerData["length"].asNumber().max()
			+ self.ribosomeElongationRate.asNumber(units.aa / units.s)
			)

		self.translationSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.translationSequences.fill(PAD_VALUE)

		aaIDs_singleLetter = self.aaIDs_singleLetter[:]

		aaMapping = {aa:i for i, aa in enumerate(aaIDs_singleLetter)}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.translationSequences[i, j] = aaMapping[letter]

		self.translationMonomerWeights = (
			(
				self.getMass(self.aaIDs)
				- self.getMass(["H2O[c]"])
				)
			/ self.nAvogadro
			).asNumber(units.fg)

		self.translationEndWeight = (self.getMass(["H2O[c]"]) / self.nAvogadro).asNumber(units.fg)

	def _buildConstants(self):
		self.constants = self._constantData
		self.__dict__.update(self.constants)


	def _buildParameters(self):
		self.__dict__.update(self._parameterData)


	def _buildMoleculeGroups(self):
		moleculeGroups = {
			'ntpIds'			:	["ATP[c]","CTP[c]","GTP[c]","UTP[c]"],
			'dNtpIds'			:	["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"],
			'dNmpIds'			:	["DAMP[c]", "DCMP[c]", "DGMP[c]", "DTMP[c]"],
			'dNmpNuclearIds'	:	["DAMP[n]", "DCMP[n]", "DGMP[n]", "DTMP[n]"],
			'rnapIds'			:	["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"],
			'polymerizedAA_IDs'	:	self._polymerizedAA_IDs, # TODO: end weight
			'polymerizedNT_IDs'	:	self._polymerizedNT_IDs, # TODO: end weight
			'polymerizedDNT_IDs':	self._polymerizedDNT_IDs,
		}

		self.__dict__.update(moleculeGroups)

	def _buildAllMasses(self):
		size = len(self._rnas) + len(self._proteins) + len(self._proteinComplexes) + len(self._metabolites) + len(self._polymerized)
		allMass = np.empty(size,
			dtype = [
					('id',		'a50'),
					('mass',	"f8")
					]
			)

		listMass = []
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self._rnas])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self._proteins])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self._proteinComplexes])
		listMass.extend([(x['id'],np.sum(x['mw7.2'])) for x in self._metabolites])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self._polymerized])

		allMass[:] = listMass

		field_units = {
			'id'		:	None,
			'mass'		:	units.g / units.mol,
			}

		self._allMass = UnitStructArray(allMass, field_units)


	def _buildTrnaData(self):
		self.aa_trna_groups = AA_TRNA_GROUPS
		self.aa_synthetase_groups = AA_SYNTHETASE_GROUPS

		# tRNA synthetase rates
		trna_synthetase_rates = SYNTHETASE_RATE
		## If no rate curated fill in with average non-zero rate
		# mean_rate = np.mean([x for x in trna_synthetase_rates.itervalues() if x != None])
		# for x in trna_synthetase_rates.iterkeys():
		# 	if trna_synthetase_rates[x] == None:
		# 		trna_synthetase_rates[x] = mean_rate
		self.trna_synthetase_rates = trna_synthetase_rates.values()
		# TODO: Remove here only for fitting
		self.synthetase_counts = None
		self.synthetase_variance = None
		self.initial_aa_polymerization_rate = None
		self.minimum_trna_synthetase_rates = None

## -- Utility functions -- ##
	def _checkDatabaseAccess(self, table):
		if len(table.objects.all()) <= 0:
			raise Exception, "Database Access Error: Cannot access public_{} table".format(table.__name__.lower())

	def getMass(self, ids):
		assert isinstance(ids, list) or isinstance(ids, np.ndarray)
		idx = [np.where(self._allMass['id'] == re.sub("\[[a-z]\]","", i))[0][0] for i in ids]
		return self._allMass['mass'][idx]
