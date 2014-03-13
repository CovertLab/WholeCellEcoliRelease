"""
Test KnowledgeBase.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

from __future__ import division

import unittest
import warnings

import cPickle
import os

import wholecell.utils.config
TEST_FIXTURE_DIR = wholecell.utils.config.TEST_FIXTURE_DIR

import nose.plugins.attrib as noseAttrib

class Test_KnowledgeBase(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		import wholecell.utils.knowledgebase_fixture_manager
		self.kb = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(os.path.join(TEST_FIXTURE_DIR, 'KnowledgeBase.cPickle'))

		# To load from the "raw" data, uncomment the following:
		# self.kb = wholecell.knowledgebase.knowledgebase.KnowledgeBase(dataFileName = "data/KnowledgeBase.xlsx",
		# 													 seqFileName = "data/KnowledgeBase.fna")

	def tearDown(self):
		pass