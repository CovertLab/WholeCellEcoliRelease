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
import Bio.SeqIO
import Bio.Seq
import Bio.Alphabet.IUPAC
import Bio.SeqUtils.ProtParam

import ipdb

class KnowledgeBaseValidator(object):
	""" KnowledgeBaseValidator """

	def __init__(self, knowledgeBase):
		self.kb = knowledgeBase

		# Validate datatypes
		self.validateMetabolites()
		self.validateProteins()
		self.validateRnas()
		self.validateGenes()
		self.validateReactions()


		# Other


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

		# Validate that biomassLoc is actual allowed location abbreviation
		s += self.checkAllowedLocation([x for x in self.kb.metabolites if x['biomassConc'] != 0.], 'biomassLoc')

		# Validate that equivEnzIds if they exist are actual proteins

		# Validate that MW is correct

		if len(s): raise Exception, s


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
						 'modifiedForm': [bool],
						 'modifiedForms': [list],
						 'monomer': [bool],
						 'mw': [float],
						 'name': [str],
						 'ntCount': [numpy.ndarray],
						 'rnaId': [str],
						 'seq': [str, unicode],
						 'unmodifiedForm': [None, str]}
		s += self.validateDatatype(fieldDataType, self.kb.proteins)

		# Validation that location that is valid
		s += self.checkAllowedLocation(self.kb.proteins, 'location')

		# Validate that composition has length >1 if it is a complex

		# Check that proteins have a geneId that exists
		# Check that complexes have no geneId
		# Check that modified forms have no geneId

		# Validate that proteins have rna taht exist

		# Validate that modified forms exist

		# Validate that unmodified forms exist

		# Validate the MW is correct and there

		# validate that sequence uses correct alphabet

	def validateRnas(self):
		s = ''
		# Validate datatypes
		fieldDataType = {'composition': [list],
						 'expression': [float],
						 'geneId': [str],
						 'halfLife': [float],
						 'id': [str, unicode],
						 'location': [str],
						 'modifiedForm': [bool],
						 'modifiedForms': [list],
						 'monomer': [bool],
						 'monomerId': [None, str, unicode],
						 'mw': [float],
						 'name': [str],
						 'ntCount': [numpy.ndarray],
						 'seq': [str],
						 'unmodifiedForm': [None, str]}
		s += self.validateDatatype(fieldDataType, self.kb.rnas)

		# Check that location is an allowed abbreviation
		s += self.checkAllowedLocation(self.kb.rnas, 'location')

		# Check that monomerId is aa legit protein id

		# If it is an mRNA then it should have a monomer id

		# Validate that expression is >0
		s = ''
		for rna in self.kb.rnas:
			if rna['expression'] < 0:
				s += '%s has a expression that is less than zero!\n' % rna['id']

		if abs(sum(rna['expression'] for rna in self.kb.rnas) - 1.0) > 1e-8:
			s += 'Expression does not sum to unity it sums to %s!\n' % str(exp_sum) 

		# Validate that half life is >0
		s = ''
		for rna in self.kb.rnas:
			if rna['halfLife'] < 0:
				s += '%s has a half life that is less than zero!\n' % rna['id']


		# Validate that its modified forms are actual frame ids

		# Validate that it has a MW and that it is correct

		# Validate that NT count is correct size

		# Validate the NT count sums to lenght of sequence

		# Validate that unmodified form is a legit frame id

		# Validate that sequence uses correct alphabet

		# Check that modified forms have no expression level

		if len(s): raise Exception, s

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

		# Validate that its rnaId is a legit one
		self.checkFrameId(self.kb.genes, 'rnaId', self.kb.rnas)

		# Validate sequence alphabet

		# Validate length against actual length of sequence
		for gene in self.kb.genes:
			if gene['length'] != len(gene['seq']):
				s += 'Gene %s has an incorrect length or sequence!\n' % gene['id']

		# Validate type is either mRNA, rRNa, tRNA, miscRNA, etc.
		for gene in self.kb.genes:
			pass

		# Check that sequence is coding strand for mRNAs
		for gene in [x for x in self.kb.genes if x['type'] == 'mRNA']:
			if gene['seq'][:3] not in ['ATG', 'GTG', 'TTG', 'ATT', 'CTG']:
				print 'Warning: Gene %s has sequence that is not the coding strand! Sequence starts with %s. May not be a coding strand mRNA.' % (gene['id'], gene['seq'][:3])

		if len(s): raise Exception, s

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

		if len(s): raise Exception, s

	def validateDatatype(self, fieldDataType, objList):
		s = ''
		for obj in objList:
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

		if len(s): return s
		else: return ''

	def checkAllowedLocation(self, listToCheck, fieldName):
		allowedLocations = [x['abbrev'] for x in self.kb.compartments]
		s = ''
		for obj in listToCheck:
			if obj[fieldName] not in allowedLocations:
				s += '%s has an invalid location abbreviation %s!\n' % (obj['id'], obj[fieldName])
		if len(s): return s
		else: return ''

	def checkFrameId(self, listToCheck, fieldName, listToCheckAgainst):
		validFrameIds = [x['id'] for x in listToCheckAgainst]
		s = ''
		for obj in listToCheck:
			if obj[fieldName] not in validFrameIds:
				s += '%s has an invalid frameId in field %s with value %s!\n' % (obj['id'], fieldName, obj[fieldName])
		if len(s): raise Exception, s


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


