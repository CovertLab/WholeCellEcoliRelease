#!/usr/bin/env python

"""
Test_BulkMolecules_partition.py

@author: Nick Ruggero, John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2013
"""

import unittest
import nose.plugins.attrib as noseAttrib
import cPickle
import os

import numpy as np
import wholecell.states.molecule_counts as wcMoleculeCounts

class Test_BulkMolecules_partition(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.countsBulkRequested = np.zeros((7, 1, 3))
		self.countsBulkRequested[..., 0] = np.array([3., 0., 0., 0., 0., 0., 0.], ndmin = 2).T
		self.countsBulkRequested[..., 1] = np.array([5., 3., 2., 2., 0., 0., 0.], ndmin = 2).T
		self.countsBulkRequested[..., 2] = np.array([20., 1., 3., 1., 2., 0., 2.], ndmin = 2).T

		self.countsBulk = np.zeros((7,1))
		self.countsBulk[:,0] = np.array([10.,2.,5.,7.,20.,3.,7.]).T
		self.countsBulkPartitioned = np.zeros((7,1,3), dtype = float)

		self.countsBulkPartitioned = np.zeros((7,1,3), dtype = float)

	def tearDown(self):
		pass

	@noseAttrib.attr('smalltest')
	def test_relativeAllocation(self):
		'''
		Tests that relative allocation works. All partitions are of lower priority.
		'''
		# isRequestAbsolute = [False, False, True, True] = 1 x # of partitoins
		# countsBulkRequested = # species x # compartments x # partitions
		# countsBulk = # speces x # compartments
		# countsBulkPartitioned = # species x # compartments x # partitions

		isRequestAbsolute = np.array([False, False, False])
		countsBulkRequested = self.countsBulkRequested
		countsBulk = self.countsBulk
		countsBulkPartitioned = self.countsBulkPartitioned

		wcMoleculeCounts.calculatePartition(isRequestAbsolute, countsBulkRequested, countsBulk, countsBulkPartitioned)

		countsBulkPartitioned_test = np.zeros((7,1,3), dtype = float)
		countsBulkPartitioned_test[...,0] = np.array([1., 0., 0., 0., 0., 0., 0.], ndmin = 2).T
		countsBulkPartitioned_test[...,1] = np.array([1., 1., 2., 4., 0., 0., 0.], ndmin = 2).T
		countsBulkPartitioned_test[...,2] = np.array([7., 0., 3., 2., 20., 0., 7.], ndmin = 2).T

		self.assertEqual(countsBulkPartitioned.tolist(), countsBulkPartitioned_test.tolist())

	@noseAttrib.attr('smalltest')
	def test_absoluteAllocation(self):
		'''
		Tests that relative allocation works with one higher priority partition.
		'''
		isRequestAbsolute = np.array([False, True, False])
		countsBulkRequested = self.countsBulkRequested
		countsBulk = self.countsBulk
		countsBulkPartitioned = self.countsBulkPartitioned

		wcMoleculeCounts.calculatePartition(isRequestAbsolute, countsBulkRequested, countsBulk, countsBulkPartitioned)

		countsBulkPartitioned_test = np.zeros((7,1,3), dtype = float)
		countsBulkPartitioned_test[...,0] = np.array([0., 0., 0., 0., 0., 0., 0.], ndmin = 2).T
		countsBulkPartitioned_test[...,1] = np.array([5., 2., 2., 2., 0., 0., 0.], ndmin = 2).T
		countsBulkPartitioned_test[...,2] = np.array([4., 0., 3., 5., 20., 0., 7.], ndmin = 2).T

		self.assertEqual(countsBulkPartitioned.tolist(), countsBulkPartitioned_test.tolist())

	@noseAttrib.attr('smalltest')
	def test_absoluteAllocation_withConflict(self):
		'''
		Tests that if two partitions are of higher priority and conflict that partitioning still works.
		'''
		isRequestAbsolute = np.array([False, True, True])
		countsBulkRequested = self.countsBulkRequested
		countsBulk = self.countsBulk
		countsBulkPartitioned = self.countsBulkPartitioned

		wcMoleculeCounts.calculatePartition(isRequestAbsolute, countsBulkRequested, countsBulk, countsBulkPartitioned)

		countsBulkPartitioned_test = np.zeros((7,1,3), dtype = float)
		countsBulkPartitioned_test[...,0] = np.array([0., 0., 0., 0., 0., 0., 0.], ndmin = 2).T
		countsBulkPartitioned_test[...,1] = np.array([2., 1., 2., 2., 0., 0., 0.], ndmin = 2).T
		countsBulkPartitioned_test[...,2] = np.array([8., 0., 3., 1., 2., 0., 2.], ndmin = 2).T

		self.assertEqual(countsBulkPartitioned.tolist(), countsBulkPartitioned_test.tolist())