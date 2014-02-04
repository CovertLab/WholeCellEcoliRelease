"""
Test KnowledgeBase.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

import unittest
import warnings

import numpy
import cPickle
import os
import wholecell.kb.KnowledgeBase

class Test_KnowledgeBase(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.kb = cPickle.load(open(os.path.join("data", "fixtures", "KnowledgeBase.cPickle"), "r"))

		# To load from the "raw" data, uncomment the following:
		# self.kb = wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileName = "data/KnowledgeBase.xlsx",
		# 													 seqFileName = "data/KnowledgeBase.fna")

	def tearDown(self):
		pass

	def test_construction(self):
		wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileName = "data/KnowledgeBase.xlsx",
													seqFileName = "data/KnowledgeBase.fna")