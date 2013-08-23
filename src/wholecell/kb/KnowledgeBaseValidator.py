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

		# Validate datatypes
		self.validateProteins()


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

	def validateProteins(self):
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
		self.validateDatatype(fieldDataType, self.kb.proteins)

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

		if len(s):
			raise Exception, s

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

	# def validateMetabolites(self):
	# 	pass


