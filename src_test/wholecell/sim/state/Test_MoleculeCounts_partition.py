#!/usr/bin/env python

"""
Test_MoleculeCounts_partition.py

@author: Nick Ruggero, John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2013
"""

import unittest
import nose.plugins.attrib as noseAttrib
import cPickle
import os

import numpy
import wholecell.sim.state.MoleculeCounts as wcMoleculeCounts

class Test_MoleculeCounts_partition(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		self.kb = cPickle.load(open(os.path.join("data","fixtures","KnowledgeBase.cPickle"), "r"))

		self.mc = self.sim.states['MoleculeCounts']
		
		self.metIds = [mol['id'] + '[c]' for mol in self.kb.metabolites if 'met' in mol['id']]

		self.metCounts = [10.,2.,5.,7.,20.,3.,7.]

	def tearDown(self):
		pass

	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
	def test_relativeAllocation(self):
		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.mc.prepartition()
		self.mc.partition()

		self.assertEqual(self.partition1.countsBulk().flatten().tolist(), [1., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.partition2.countsBulk().flatten().tolist(), [1., 1., 2., 4., 0., 0., 0.])
		self.assertEqual(self.partition3.countsBulk().flatten().tolist(), [7., 0., 3., 2., 20., 0., 7.])

	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
	def test_absoluteAllocation(self):
		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.partition2.isReqAbs = True # this is a hack

		self.mc.prepartition()
		self.mc.partition()
		
		self.assertEqual(self.partition1.countsBulk().flatten().tolist(), [0., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.partition2.countsBulk().flatten().tolist(), [5., 2., 2., 2., 0., 0., 0.])
		self.assertEqual(self.partition3.countsBulk().flatten().tolist(), [4., 0., 3., 5., 20., 0., 7.])

	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
	def test_absoluteAllocation_withConflict(self):
		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.partition2.isReqAbs = True # this is a hack
		self.partition3.isReqAbs = True # this is a hack

		self.mc.prepartition()
		self.mc.partition()
		
		self.assertEqual(self.partition1.countsBulk().flatten().tolist(), [0., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.partition2.countsBulk().flatten().tolist(), [2., 1., 2., 2., 0., 0., 0.])
		self.assertEqual(self.partition3.countsBulk().flatten().tolist(), [8., 0., 3., 1., 2., 0., 2.])