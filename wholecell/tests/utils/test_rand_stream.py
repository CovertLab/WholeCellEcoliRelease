#!/usr/bin/env python

"""
Test randStream.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

from __future__ import division

import unittest
import warnings

import numpy as np
import wholecell.utils.rand_stream

import nose.plugins.attrib as noseAttrib

class Test_randStream(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		cls.randStream = wholecell.utils.rand_stream.RandStream(seed = 1)

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.randStream.reset(1)

	def tearDown(self):
		pass

	@noseAttrib.attr('smalltest')
	def test_randsample(self):
		r = self.randStream

		with self.assertRaises(AttributeError) as context:
			r.randw(np.array([[0, 0], [0, 0]]))
		self.assertEqual(context.exception.message, "'RandStream' object has no attribute 'randw'")

		# NOTE: In Matlab, no exception is thrown here (the returned array is merely empty)
		with self.assertRaises(ValueError) as context:
			self.assertEqual(np.zeros([0, 1]), r.randsample(2, 1, True, [0, 0]))
		self.assertEqual(context.exception.message, "weights sum to zero - cannot scale")

		with self.assertRaises(ValueError) as context:
			r.randsample(2, 1, True, [-1, 0])
		self.assertEqual(context.exception.message, "probabilities are not non-negative")

		with self.assertRaises(ValueError) as context:
			r.randsample(2, 1, True, [0, np.Inf])
		self.assertEqual(context.exception.message, "expect finite weights")

		with warnings.catch_warnings(record = True) as w:
			with self.assertRaises(ValueError) as context:
				r.randsample(2, 1, True, [0, np.NaN])
			self.assertEqual(context.exception.message, "probabilities do not sum to 1")

		with self.assertRaises(ValueError) as context:
			r.randsample(2, -1, True, [0, 1])
		self.assertEqual(context.exception.message, "negative dimensions are not allowed")

		# NOTE: In Matlab, there is a test for integer-valued k
		# However, np will just floor a float-valued k and use it

		with self.assertRaises(OverflowError) as context:
			r.randsample(2, np.Inf, True, [0, 1])
		self.assertEqual(context.exception.message, "cannot convert float infinity to integer")

		with self.assertRaises(ValueError) as context:
			r.randsample(2, np.NaN, True, [0, 1])
		self.assertEqual(context.exception.message, "cannot convert float NaN to integer")

		with self.assertRaises(Exception) as context:
			r.randsample(2, [1, 1], True, [0, 1])
		self.assertEqual(context.exception.message, "Expect k to be a scalar")

		# Sampling
		self.assertEqual(1, r.randsample(7, 1, True, [0, 1, 1, 1, 1, 0, 2]).size)
		self.assertEqual(0, r.randsample(7, 0, True, [0, 1, 1, 1, 1, 0, 2]).size)
		self.assertEqual(10, r.randsample(7, 10, True, [0, 1, 1, 1, 1, 0, 2]).size)

		with self.assertRaises(ValueError) as context:
			r.randsample(7, 10, False, [0, 1, 1, 1, 1, 0, 2])
		self.assertEqual(context.exception.message, "Cannot take a larger sample than population when 'replace=False'")

		self.assertTrue(np.array_equal(np.array([6, 2, 1, 3, 4]), r.randsample(7, 5, False, [0, 1, 1, 1, 1, 0, 2])))
		self.assertTrue(np.array_equal(np.array([1, 2, 6, 4, 3]), r.randsample(7, 5, False, np.transpose([0, 1, 1, 1, 1, 0, 2]))))

		# Fairness
		counts = np.zeros([2, 1])
		for i in xrange(100):
			idx = r.randsample(2, 1, True, [1, 1])
			counts[idx] += 1
		self.assertTrue(np.max(counts) - np.min(counts) < 0.12 * np.max(counts))

	@noseAttrib.attr('smalltest')
	def test_randomlySelectRows(self):
		r = self.randStream

		mat = np.reshape(np.arange(10), [5, 2])
		self.assertTrue(np.array_equal(mat, r.randomlySelectRows(mat, 1)[0]))

		mat = np.reshape(np.arange(4), [2, 2])
		counts = np.zeros([2, 1])
		for i in xrange(1000):
			selMat, selIdx = r.randomlySelectRows(mat, 0.5)
			self.assertTrue(np.array_equal(selMat, mat[selIdx, :]), "Randomly selected rows should correspond to randomly selected indices")
			self.assertTrue(np.array_equal(selIdx, np.unique(selIdx)), "Random selection should not be done with replacement")
			counts[selIdx] += 1
		self.assertTrue(np.max(counts) - np.min(counts) < 0.11 * np.max(counts), "Random selection should be unbiased")

	@noseAttrib.attr('smalltest')
	def test_randomlySelectNRows(self):
		r = self.randStream

		mat = np.reshape(np.arange(10), [5, 2])
		self.assertTrue(np.array_equal(mat, r.randomlySelectRows(mat, 1)[0]))
		self.assertEqual(5, r.randomlySelectNRows(mat, 5)[0].shape[0])

		mat = np.reshape(np.arange(4), [2, 2])
		counts = np.zeros(2)
		for i in xrange(1000):
			selMat, selIdx = r.randomlySelectNRows(mat, 1)
			self.assertTrue(np.array_equal(selMat, mat[selIdx, :]), "Randomly selected rows should correspond to randomly selected indices")
			self.assertTrue(np.array_equal(selIdx, np.unique(selIdx)), "Random selection should not be done with replacement")
			counts[selIdx] += 1
		self.assertTrue(np.max(counts) - np.min(counts) < 0.11 * np.max(counts), "Random selection should be unbiased")

	@noseAttrib.attr('smalltest')
	def test_randCounts(self):
		r = self.randStream

		with self.assertRaises(Exception) as context:
			r.randCounts(2, -1)
		self.assertEqual(context.exception.message, "N must be positive.")

		self.assertEqual(0, r.randCounts(2, 0))
		self.assertEqual(1, r.randCounts(2, 1))
		self.assertEqual(2, r.randCounts(2, 2))

		self.assertTrue(np.array_equal([2, 2, 3], r.randCounts([2, 2, 3], 7)))
		with self.assertRaises(Exception) as context:
			r.randCounts([2, 2, 3], 8)
		self.assertEqual(context.exception.message, "N must be at most the total available counts.")	

		self.assertEqual(2, np.sum(r.randCounts([2, 2, 3], 2)))

		self.assertTrue(np.all([2, 2, 3] >= r.randCounts([2, 2, 3], 6)))

		counts = np.zeros(2)
		for i in xrange(1000):

			tmp = r.randCounts([2, 2], 1)
			self.assertEqual(1, np.sum(tmp))
			counts += tmp

		self.assertTrue(np.max(counts) - np.min(counts) < 0.10 * np.max(counts), "Random selection should be unbiased")

	@noseAttrib.attr('smalltest')
	def test_getters(self):
		r = self.randStream
		self.assertEqual('MT19937', r.type)
		self.assertEqual(5, len(r.state))

	@noseAttrib.attr('smalltest')
	def test_setters(self):
		r = self.randStream
		state = r.state
		r.randi(5)
		self.assertFalse(np.array_equal(state[1], r.state[1]))
		r.state = state
		self.assertTrue(np.array_equal(state[1], r.state[1]))
		r.seed = 1234
		self.assertFalse(np.array_equal(state[1], r.state[1]))