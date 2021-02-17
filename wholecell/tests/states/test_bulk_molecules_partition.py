"""
Test_BulkMolecules_partition.py
"""

from __future__ import absolute_import, division, print_function

import unittest

import numpy as np
import wholecell.states.bulk_molecules as wcBulkMolecules


class Test_BulkMolecules_partition(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.countsBulkRequested = np.zeros((7, 3))
		self.countsBulkRequested[..., 0] = np.array([3., 0., 0., 0., 0., 0., 0.]).T
		self.countsBulkRequested[..., 1] = np.array([5., 3., 2., 2., 0., 0., 0.]).T
		self.countsBulkRequested[..., 2] = np.array([20., 1., 3., 1., 2., 0., 2.]).T

		self.countsBulk = np.array([10.,2.,5.,7.,20.,3.,7.])
		self.countsBulkPartitioned = np.zeros((7, 3), dtype = float)

	def tearDown(self):
		pass

	def test_relativeAllocation(self):
		'''
		Tests that relative allocation works. All partitions are of lower priority.
		'''
		# processPriorities = [False, False, True, True] = 1 x # of partitoins
		# countsBulkRequested = # species x # compartments x # partitions
		# countsBulk = # speces x # compartments
		# countsBulkPartitioned = # species x # compartments x # partitions

		random_state = np.random.RandomState(0)
		processPriorities = np.array([0, 0, 0])
		countsBulkRequested = self.countsBulkRequested
		countsBulk = self.countsBulk

		countsBulkPartitioned = wcBulkMolecules.calculatePartition(
			processPriorities, countsBulkRequested, countsBulk, random_state)

		countsBulkPartitioned_test = np.zeros((7, 3), dtype = float)
		countsBulkPartitioned_test[...,0] = np.array([1., 0., 0., 0., 0., 0., 0.]).T
		countsBulkPartitioned_test[...,1] = np.array([2., 1., 2., 2., 0., 0., 0.]).T
		countsBulkPartitioned_test[...,2] = np.array([7., 1., 3., 1., 2., 0., 2.]).T

		self.assertEqual(countsBulkPartitioned.tolist(), countsBulkPartitioned_test.tolist())

	def test_absoluteAllocation(self):
		'''
		Tests that relative allocation works with one higher priority partition.
		'''
		random_state = np.random.RandomState(0)
		processPriorities = np.array([0, 10, 0])
		countsBulkRequested = self.countsBulkRequested
		countsBulk = self.countsBulk

		countsBulkPartitioned = wcBulkMolecules.calculatePartition(
			processPriorities, countsBulkRequested, countsBulk, random_state)

		countsBulkPartitioned_test = np.zeros((7, 3), dtype = float)
		countsBulkPartitioned_test[...,0] = np.array([1., 0., 0., 0., 0., 0., 0.]).T
		countsBulkPartitioned_test[...,1] = np.array([5., 2., 2., 2., 0., 0., 0.]).T
		countsBulkPartitioned_test[...,2] = np.array([4., 0., 3., 1., 2., 0., 2.]).T

		self.assertEqual(countsBulkPartitioned.tolist(), countsBulkPartitioned_test.tolist())

	def test_absoluteAllocation_withConflict(self):
		'''
		Tests that if two partitions are of higher priority and conflict that partitioning still works.
		'''
		random_state = np.random.RandomState(0)
		processPriorities = np.array([0, 10, 10])
		countsBulkRequested = self.countsBulkRequested
		countsBulk = self.countsBulk

		countsBulkPartitioned = wcBulkMolecules.calculatePartition(
			processPriorities, countsBulkRequested, countsBulk, random_state)

		countsBulkPartitioned_test = np.zeros((7, 3), dtype = float)
		countsBulkPartitioned_test[...,0] = np.array([0., 0., 0., 0., 0., 0., 0.]).T
		countsBulkPartitioned_test[...,1] = np.array([2., 1., 2., 2., 0., 0., 0.]).T
		countsBulkPartitioned_test[...,2] = np.array([8., 1., 3., 1., 2., 0., 2.]).T

		self.assertEqual(countsBulkPartitioned.tolist(), countsBulkPartitioned_test.tolist())

	def test_high_counts(self):
		"""
		Tests specific counts case where rounding caused a ValueError in previous versions.
		"""
		random_state = np.random.RandomState(0)
		processPriorities = np.array([0, 0, 0, 0, 0, 0, 0])
		countsBulkRequested = np.array([[0, 14768, 0, 609520, 0, 0, 0]])
		countsBulk = np.array([604314])

		wcBulkMolecules.calculatePartition(
			processPriorities, countsBulkRequested, countsBulk, random_state)

	def test_excess_remainder(self):
		"""
		Tests that a remainder greater than 1 will be handled without an error.
		"""
		random_state = np.random.RandomState(0)
		processPriorities = np.array([0, 0, 0])
		countsBulkRequested = np.array([[7, 2, 1]])
		countsBulk = np.array([4])

		wcBulkMolecules.calculatePartition(
			processPriorities, countsBulkRequested, countsBulk, random_state)
