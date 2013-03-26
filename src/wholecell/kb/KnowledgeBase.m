#!/usr/bin/env python

"""
KnowledgeBase

Whole-cell knowledge base

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

import os.path
import openpyxl as xl


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
		pass

	def loadGenome(self):
		pass

	def loadGenes(self):
		pass

	def loadComplexes(self):
		pass

	def loadReactions(self):
		pass