#!/usr/bin/env python

"""
Test_BulkMolecules.py

@author: Nick Ruggero, John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2013
"""

from __future__ import division

import unittest
import nose.plugins.attrib as noseAttrib
import cPickle
import os

import numpy as np
import wholecell.states.molecule_counts as wcMoleculeCounts

# class Test_BulkMolecules(unittest.TestCase):

# 	@classmethod
# 	def setUpClass(cls):
# 		pass

# 	@classmethod
# 	def tearDownClass(cls):
# 		pass

# 	def setUp(self):
# 		self.sim = cPickle.load(open(os.path.join("fixtures", "Simulation.cPickle"), "r"))
# 		self.kb = cPickle.load(open(os.path.join("fixtures","KnowledgeBase.cPickle"), "r"))

# 		self.bulkMolecules = self.sim.states['BulkMolecules']

# 		self.rand_counts = np.random.randint(0, 2000, size = self.bulkMolecules.countsBulk().shape)
# 		self.bulkMolecules.countsBulkIs(self.rand_counts)

# 	def tearDown(self):
# 		pass

# 	@noseAttrib.attr('smalltest')
# 	def test_CountsBulkViews(self):
# 		view1 = self.bulkMolecules.countsBulkViewNew()
# 		view2 = self.bulkMolecules.countsBulkViewNew([self.bulkMolecules._molIDs[0]
# 					 + '[' + self.bulkMolecules._compartments[0] + ']',])

# 		# Test accessing
# 		self.assertEqual(view1.countsBulk().tolist(),
# 			self.rand_counts.tolist())
		
# 		# Test modification within a view
# 		view1.countsBulkIs(view1.countsBulk() + 1)
# 		self.assertEqual((view1.countsBulk() - 1).tolist(),
# 			self.rand_counts.tolist())

# 		# Test modification across views
# 		self.assertEqual(view2.countsBulk() - 1, self.rand_counts[0,0])

# 	@noseAttrib.attr('smalltest')
# 	def test_CountsBulk(self):
# 		self.assertEqual(
# 			self.bulkMolecules.countsBulk().tolist(),
# 			self.rand_counts.tolist()
# 			)

# 	@noseAttrib.attr('smalltest')
# 	def test_PartitionIndexing(self):
# 		metabolite = 'ATP'
# 		compartment = 'c'

# 		partition_of_interest = self.bulkMolecules.partitions['Transcription']
# 		molIdx, cmpIdx = partition_of_interest._getIndex(
# 			'{}[{}]'.format(metabolite, compartment)
# 			)[1:]

# 		# Insure that _getIndex points to the right metabolite
# 		self.assertEqual(metabolite, partition_of_interest._molIDs[molIdx])
# 		self.assertEqual('merged', partition_of_interest._compartments[cmpIdx])

# 		# TODO: add more tests, more realistic tests cases
