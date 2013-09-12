#!/usr/bin/env python

"""
KnowledgeBaseValidator

Validates Whole-cell knowledge base

@author: Nicholas Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/22/2013
"""

import os.path
import openpyxl as xl
import csv
import json
import numpy
import itertools
import inspect

import ipdb

class KnowledgeBaseValidator(object):
	""" KnowledgeBaseValidator """

	def __init__(self, knowledgeBase):
		self.kb = knowledgeBase

		s = ''
		# Validate datatypes
		s += self.validateMetabolites()
		s += self.validateProteins()
		s += self.validateRnas()
		s += self.validateGenes()
		s += self.validateReactions()


		# Other

		if len(s):
			raise Exception, s

	def validateMetabolites(self):
		s = ''
		# Validate datatypes
		fieldDataType = {'biomassConc': [float],
						 'biomassLoc': [str, None],
						 'charge7.2': [int],
						 'comments': [str],
						 'equivEnzIds': [None, list],
						 'fakeMet': [bool],
						 'formula7.2': [str],
						 'formulaNeutral': [str],
						 'id': [str, unicode],
						 'maxExchange': [float],
						 'mediaConc': [float],
						 'mw7.2': [float],
						 'name': [str]}
		s += self.validateDatatype(fieldDataType, self.kb.metabolites)

		metDict = dict([(x['id'],x) for x in self.kb.metabolites])

		# Check that biomassConc was >= 0
		for met in self.kb.metabolites:
			if met['biomassConc'] < 0:
				s += 'Metabolite %s has a biomassConc <0!\n' % met['id']

		# Validate that biomassLoc is actual allowed location abbreviation
		s += self.checkAllowedLocation([x for x in self.kb.metabolites if x['biomassConc'] != 0.], 'biomassLoc')

		# Validate that equivEnzIds if they exist are actual proteins and their location is valid
		validProteinIds = [x['id'] for x in self.kb.proteins]
		validRnas = [x['id'] for x in self.kb.rnas]
		allowedLocations = [x['abbrev'] for x in self.kb.compartments]
		for met in self.kb.metabolites:
			for equivEnz in met['equivEnzIds']:
				if not equivEnz['id'] in validProteinIds and not equivEnz['id'] in validRnas:
					s += 'Fake metabolite %s has an enzyme id %s that is not valid!\n' % (met['id'], equivEnz['id'])
				if not equivEnz['location'] in allowedLocations:
					s += 'Fake metabolite %s has an enzyme location %s that is not valid!\n' % (met['id'], equivEnz['id'])
			if met['fakeMet'] == False and len(met['equivEnzIds']):
				s += 'Real metabolite %s has psudo-metabolite equivalents!\n' % met['id']
			if met['fakeMet'] == True and not len(met['equivEnzIds']):
				s += 'Psudo-metabolite %s has no equivalent protein ids!\n' % met['id']

		# Validate that MW is correct for key species and assume the rest is calculated correctly
		testDict = {'H'		: 1.0079,
					'H2O'	: 18.0148,
					'CYS-L'	: 121.1533,
					'SELNP'	: 159.9468}

		for metId in testDict.iterkeys():
			if abs(metDict[metId]['mw7.2'] - testDict[metId]) > 1e-8:
				s += 'Molecular weights are not calculated correctly for metabolites! Test case %s\n' % metId

		return s

	def validateProteins(self):
		s = ''
		# Validate datatypes
		fieldDataType = {'aaCount': [numpy.ndarray],
						 'comments': [str],
						 'composition': [list],
						 'formationProcess': [str],
						 'geneId': [str],
						 'id': [str, unicode],
						 'location': [str],
						 'modifiedForms': [list],
						 'mw': [float],
						 'name': [str],
						 'ntCount': [numpy.ndarray],
						 'rnaId': [str],
						 'seq': [str, unicode],
						 'unmodifiedForm': [None, str]}
		s += self.validateDatatype(fieldDataType, self.kb.proteins)

		# Validation that location that is valid
		s += self.checkAllowedLocation(self.kb.proteins, 'location')

		## Validate protein monomers that are unmodified
		validGeneIds = [x['id'] for x in self.kb.genes]
		validRnaIds = [x['id'] for x in self.kb.rnas]
		validProteinIds = [x['id'] for x in self.kb.proteins]
		rnaDict = dict([(x['id'], x) for x in self.kb.rnas])
		geneDict = dict([(x['id'], x) for x in self.kb.genes])
		protDict = dict([(x['id'], x) for x in self.kb.proteins])
		for prot in [x for x in self.kb.proteins if not len(x['composition']) and x['unmodifiedForm'] == None]:
			# Check that proteins have a geneId that exists
			if not prot['geneId'] in validGeneIds:
				s += 'Protein %s has an invalid gene id %s!\n' % (prot['id'], prot['geneId'])
			# Validate that proteins have rna that exist and that rna points to this protein and that rna points to this gene
			if not prot['rnaId'] in validRnaIds:
				s += 'Protein %s has an invalid rna id %s!\n' % (prot['id'], prot['rnaId'])
			if rnaDict[prot['rnaId']]['monomerId'] != prot['id']:
				s += 'Protein %s has an rna %s that does not point back to the protein!\n' % (prot['id'], prot['rnaId'])
			if rnaDict[prot['rnaId']]['geneId'] != prot['geneId']:
				s += 'Protein %s has an rna %s that does not point back to the same geneId!\n' % (prot['id'], prot['rnaId'])
			# Validate that modified forms exist and that unmodified form points back to this protein monomer
			for modForm in prot['modifiedForms']:
				if not modForm in validProteinIds:
					s += 'Protein %s has an invalid modified form %s!\n' % (prot['id'], modForm)
				if protDict[modForm]['unmodifiedForm'] != prot['id']:
					s += 'Protein %s has a modified form %s that does not point back to the protein!\n' % (prot['id'], modForm)
			# Validate that sequence uses correct alphabet
			proteinAlphabet = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "U", "S", "T", "W", "Y", "V"]
			letter_sum = 0
			for letter in proteinAlphabet:
				letter_sum += prot['seq'].count(letter)
			if letter_sum != len(prot['seq']):
				s += '%s has an invalid character in its sequence!\n' % prot['id']

			# Validate that sequence has the same sum as ntCount
			# Validate the ntCount's sum is equal to sequence length

		## Validate protein monomers that are modified
			# Check that modified protein forms have no geneId
			# Validate that unmodified forms exist


		## Validate complexes
		# Check that complexes have no geneId
		# Check that modified form exists and points back to this unmodified complex
		# Check that it has no sequence

		# Check that modified complexes have no geneId
		# Validate that unmodified forms exist and points to this modified complex
		# Validate the MW is correct and there




		# Validate protein monomers
		# Check that proteins have a geneId that exists
		# Validate that proteins have rna that exist and that rna points to this protein
		# Validate that modified forms exist and that unmodified form points back to this protein monomer
		# Validate that sequence uses correct alphabet
		# Validate that sequence has the same sum as ntCount
		# Validate the ntCount's sum is equal to sequence length

		# Check that modified protein forms have no geneId
		# Validate that unmodified forms exist

		# Validate the MW is correct and there





		return s

	def validateRnas(self):
		s = ''
		# Validate datatypes
		fieldDataType = {'composition': [list],
						 'expression': [float],
						 'geneId': [str],
						 'halfLife': [float],
						 'id': [str, unicode],
						 'location': [str],
						 'modifiedForms': [list],
						 'monomerId': [None, str, unicode],
						 'mw': [float],
						 'name': [str],
						 'ntCount': [numpy.ndarray],
						 'seq': [str],
						 'unmodifiedForm': [None, str]}
		s += self.validateDatatype(fieldDataType, self.kb.rnas)

		proteinDict = dict([(x['id'], x) for x in self.kb.proteins])
		geneDict = dict([(x['id'], x) for x in self.kb.genes])
		rnaDict = dict([(x['id'], x) for x in self.kb.rnas])

		# Validate that all RNAs have no composition
		if len([x for x in self.kb.rnas if x['composition'] != []]):
			s += 'An RNA has a composition!\n'

		# Validate that expression is >0
		for rna in self.kb.rnas:
			if rna['expression'] < 0:
				s += '%s has a expression that is less than zero!\n' % rna['id']

		if abs(sum(rna['expression'] for rna in self.kb.rnas) - 1.0) > 1e-8:
			s += 'Expression does not sum to unity it sums to %s!\n' % str(exp_sum) 

		# Validate that geneId is a valid id and that the gene points to the rna
		validGeneIds = [x['id'] for x in self.kb.genes]
		for rna in self.kb.rnas:
			if not rna['geneId'] in validGeneIds:
				s += 'RNA %s has invalid gene id %s!\n' % (rna['id'], rna['geneId'])

			if geneDict[rna['geneId']]['rnaId'] != rna['id']:
				s += 'RNA %s has gene %s that has incorrect or invalid rna pointer!\n' % (rna['id'], rna['geneId'])

		# Validate that halfLife is >0
		s = ''
		for rna in self.kb.rnas:
			if rna['halfLife'] < 0:
				s += '%s has a half life that is less than zero!\n' % rna['id']

		# Check that location is an allowed abbreviation
		s += self.checkAllowedLocation(self.kb.rnas, 'location')

		## Check modifiedForms properties
		for modRna in [x for x in self.kb.rnas if x['unmodifiedForm'] != None]:
			# Validate that unmodified form is a legit frame id
			rnaFrameIds = [x['id'] for x in self.kb.rnas]
			if not modRna['unmodifiedForm'] in rnaFrameIds:
				s += 'Modified RNA %s has invalid frame id for unmodifiedForm %s!\n' % (modRna['id'], modRna['unmodifiedForm'])

			# Check that modified forms have no expression level
			if modRna['expression'] != 0.:
				s += 'Modified RNA %s has a non-zero expression level!\n' % modRna['id']

			# Check that unmodified form has this modified RNA as its modified form
			unmodifiedRna = rnaDict[modRna['unmodifiedForm']]
			foundModForm = False
			for possibleModForm in unmodifiedRna['modifiedForms']:
				if possibleModForm == modRna['id']:
					foundModForm = True
			if not foundModForm:
				s += 'Unmodified RNA %s has incorrect or invalid frame id %s for possible modified form %s!\n' % (unmodifiedRna['id'], modRna['id'], unmodifiedRna['modifiedForms'])

		## Check monomerId properties
		mRNAs = [y for y in self.kb.genes if y['type'] == 'mRNA']
		proteins = [x['id'] for x in self.kb.proteins]
		for mRNA in [x for x in self.kb.rnas if x['id'] in mRNAs]:
 			# If it is an mRNA then it should have a monomer id
 			if mRNA['monomerId'] == None:
 				s += 'mRNA %s has no monomerId!\n' % mRNA['id']
			# Check monomer id is a legit protein id
			if mRNA['monomerId'] not in proteins:
 				s += 'mRNA %s has invalid monomerId. It is not in proteins!\n' % mRNA['id']

		# Validate that it has a MW and that it is correct
		for rna in self.kb.rnas:
			compare_mw = self.calculateRnaMW(rna['seq'])
			if not abs(compare_mw - rna['mw']) < 1e-6:
				s += 'RNA %s has a molecular weight of %s but when recalculated it should be %s!\n' % (rna['id'], str(rna['mw']), str(compare_mw))

		# Validate that ntCount is correct dimension and sums to length of sequence
		for rna in self.kb.rnas:
			if len(rna['ntCount']) != 4:
				s += 'RNA %s has ntCount that is incorrect in dimension!\n' % rna['id']
			if sum(rna['ntCount']) != len(rna['seq']):
				s += 'RNA %s has ntCount that is incorrect in length!\n' % rna['id']

		# Validate sequence alphabet
		s += self.validateAlphabet(self.kb.rnas, ['A','U','G','C'])

		return s

	def validateGenes(self):
		s = ''
		# Validate datatypes
		fieldDataType = {'coordinate': [int],
						 'direction': [str],
						 'id': [str],
						 'length': [int],
						 'name': [str],
						 'rnaId': [str],
						 'seq': [str],
						 'symbol': [str],
						 'type': [str]}
		s += self.validateDatatype(fieldDataType, self.kb.genes)

		# Validate that coordinate is in range of genome
		for gene in self.kb.genes:
			if gene['coordinate'] > len(self.kb.genomeSeq):
				s += 'Gene %s has coordinate greater than length of genome!\n' % gene['id']

		# Validate that direction is either a + or a -
		for gene in self.kb.genes:
			if gene['direction'] not in ['+', '-']:
				s += 'Gene %s has invalid direction!\n' % gene['id']

		# Validate that length is < length of genome
		for gene in self.kb.genes:
			if gene['length'] >= len(self.kb.genomeSeq):
				s += 'Gene %s has length longer than the genome!\n' % gene['id']

		# Validate length against actual length of sequence
		for gene in self.kb.genes:
			if gene['length'] != len(gene['seq']):
				s += 'Gene %s has an incorrect length or sequence!\n' % gene['id']

		# Validate that its rnaId is a legit one
		s += self.checkFrameId(self.kb.genes, 'rnaId', self.kb.rnas)

		# Validate sequence alphabet
		s += self.validateAlphabet(self.kb.genes, ['A','T','G','C'])

		# Validate type is either mRNA, rRNa, tRNA, miscRNA, etc.
		for gene in self.kb.genes:
			if gene['type'] not in ['mRNA', 'rRNA', 'tRNA', 'miscRNA']:
				s += 'Gene %s has an invalid type!\n' % gene['id']

		# Check that sequence is coding strand for mRNAs
		for gene in [x for x in self.kb.genes if x['type'] == 'mRNA']:
			if gene['seq'][:3] not in ['ATG', 'GTG', 'TTG', 'ATT', 'CTG']:
				print 'Warning: Gene %s has sequence that is not the coding strand! Sequence starts with %s. May not be a coding strand mRNA.' % (gene['id'], gene['seq'][:3])

		return s

	def validateReactions(self):
		s = ''
		fieldDataType = {'catBy': [list],
						 'dir': [int],
						 'ec': [str],
						 'id': [str, unicode],
						 'name': [str],
						 'process': [str],
						 'stoichiometry': [list]}
		s += self.validateDatatype(fieldDataType, self.kb.reactions)

		# Check that catBy actually contains list of enzymes

		# Check that direction is 1 or -1

		# Check that process name is valid

		# Check that stoichiometry datatypes are valid metabolites

		# Check that metabolties are real metabolites

		# Check for mass balance

		return s

	# def validateFrameId(self):
	# 	# Validates that all frameid's are unique
	# 	# TODO: Redo this so that it uses the id's actually in the KB
	# 	frameIds = []
	# 	fileNames = []
	# 	self.getFrameIds('genes.csv', frameIds, fileNames)
	# 	self.getFrameIds('locations.csv', frameIds, fileNames)
	# 	self.getFrameIds('metabolites.csv', frameIds, fileNames)
	# 	self.getFrameIds('promoters.csv', frameIds, fileNames)
	# 	self.getFrameIds('proteinComplexes_modified.csv', frameIds, fileNames)
	# 	self.getFrameIds('proteinComplexes.csv', frameIds, fileNames)
	# 	self.getFrameIds('proteinMonomers_modified.csv', frameIds, fileNames)
	# 	self.getFrameIds('proteinMonomers.csv', frameIds, fileNames)
	# 	self.getFrameIds('reactions.csv', frameIds, fileNames)
	# 	self.getFrameIds('rna_modified.csv', frameIds, fileNames)
	# 	self.getFrameIds('rna.csv', frameIds, fileNames)
	# 	self.getFrameIds('terminators.csv', frameIds, fileNames)
	# 	self.getFrameIds('transcriptionUnits.csv', frameIds, fileNames)

	# def getFrameIds(self, fileName, frameIds, fileNames):
	# 	fileName = self.kb.dataFileDir + os.sep + fileName
	# 	if not os.path.isfile(fileName):
	# 		raise Exception, "%s is missing" % fileName

	# 	with open(fileName, "r") as csvfile:
	# 		dr = csv.DictReader(csvfile, delimiter = "\t")
	# 		for row in dr:
	# 			if row['Frame ID'] in frameIds:
	# 				# TODO: Change back into an exception once we hear back from Ecocyc about the one error here
	# 				#raise Exception, '%s already in use as Frame ID!' % row['Frame ID']
	# 				print '%s already in use as Frame ID!' % row['Frame ID']
	# 				print 'frameid that offends is in filename %s ' % fileName
	# 				print 'frameid in conflict is in filename %s ' % fileNames[frameIds.index(row['Frame ID'])]
	# 			else:
	# 				frameIds.append(row['Frame ID'])
	# 				fileNames.append(fileName)

	def validateCenteralDogmaConnections(self):
		pass

		# Check gene-->rna-->gene product

	def validateDatatype(self, fieldDataType, objList):
		s = ''
		for obj in objList:
			# Validate field names
			for field in obj.iterkeys():
				if not fieldDataType.has_key(field):
					s += '%s has an extra field with name %s!\n' % (obj['id'], field)
			for field in fieldDataType.iterkeys():
				if not obj.has_key(field):
					s += '%s is missing field with name %s!\n' % (obj['id'], field)

			# Validate datatypes
			for fieldName in fieldDataType.iterkeys():
				hasField = False
				if obj.has_key(fieldName):
					hasField = True
				else:
					s += '%s is missing field "%s"!\n' % (obj['id'], fieldName)

				if hasField:
					isOk = False
					for allowedType in fieldDataType[fieldName]:
						if inspect.isclass(allowedType):
							if isinstance(obj[fieldName], allowedType):
								isOk = True
						elif obj[fieldName] == allowedType:
							isOk = True
					if not isOk:
						s += '%s has field "%s" that is invalid with value "%s"! Has type %s requires type %s.\n' % (obj['id'], fieldName, str(obj[fieldName]), str(type(obj[fieldName])), str(fieldDataType[fieldName]))

		return s

	def checkAllowedLocation(self, listToCheck, fieldName):
		allowedLocations = [x['abbrev'] for x in self.kb.compartments]
		s = ''
		for obj in listToCheck:
			if obj[fieldName] not in allowedLocations:
				s += '%s has an invalid location abbreviation %s!\n' % (obj['id'], obj[fieldName])
		return s

	def checkFrameId(self, listToCheck, fieldName, listToCheckAgainst):
		validFrameIds = [x['id'] for x in listToCheckAgainst]
		s = ''
		for obj in listToCheck:
			if obj[fieldName] not in validFrameIds:
				s += '%s has an invalid frameId in field %s with value %s!\n' % (obj['id'], fieldName, obj[fieldName])
		return s

	def validateAlphabet(self, listToCheck, alphabet):
		s = ''
		for obj in listToCheck:
			letter_sum = 0
			for letter in alphabet:
				letter_sum += obj['seq'].count(letter)
			if letter_sum != len(obj['seq']):
				s += '%s has an invalid character in its sequence!\n' % obj['id']
		return s

	def calculatePeptideMW(self, seq):
		# Borrowed from BioPython
		# They have the Selenocysteine (U) value commented out in IUPACData.py, so we can't use their functions
		# Fortunately they are simple functions with source code available at: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-pysrc.html
		aaWeights = {
			"A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19, "G": 75.07, "H": 155.16, "I": 131.18, "K": 146.19, "L": 131.18,
			"M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09, "T": 119.12, "U": 168.05, "V": 117.15, "W": 204.23,
			"Y": 181.19
		}
		total_weight_with_water = 0.
		for base in seq:
			total_weight_with_water += aaWeights[base]
		water = 18.02
		total_weight_without_water = total_weight_with_water - (len(seq) - 1) * water
		return total_weight_without_water

	def calculateRnaMW(self, seq):
		# Borrowed from BioPython and modified to be at pH 7.2
		rnaNtMw = { 
			"A": 345.20,
			"C": 321.18,
			"G": 361.20,
			"U": 322.17,
		}
		return sum(rnaNtMw[x] for x in seq) - (len(seq) - 1) * 17.01

	# def validateMetabolicNetwork(self):
	# 	# Validate all metabolites are used in the metabolic network
	# 	metabolites = set([m['id'] for m in self.metabolites])
	# 	metabolitesUsed = set(list(itertools.chain(*[[met['molecule'] for met in rxn['stoichiometry']] for rxn in self.reactions])))

	# 	if len(metaboltiesUsed.difference(metabolites)):
	# 		raise Exception, 'Reactions use metabolites not in KnowledgeBase\n'

	# 	# TODO: Check inverse difference. SHould be zero if you count all the subunits of complexes and stuff.
	# 	problemMet = list(metabolites.difference(metabolitesUsed))

	# 	# Validate all reactions have valid metabolites


	# 	# Validate all reactions occur in allowed compartments


