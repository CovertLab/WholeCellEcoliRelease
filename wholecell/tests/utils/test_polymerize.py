"""
Test polymerize.py

	cd wcEcoli
	pytest wholecell/tests/utils/test_polymerize.py
"""

from __future__ import absolute_import, division, print_function

from wholecell.utils.polymerize import (buildSequences, polymerize,
	computeMassIncrease, sum_monomers, sum_monomers_reference_implementation)

import numpy as np
from numpy.testing import assert_equal

import unittest
from six.moves import range

P = polymerize.PAD_VALUE


class Test_polymerize(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		pass

	def tearDown(self):
		pass

	def makeSequenceMonomers(self, sequences):
		"""
		Set self.sequences, self.nSequences, self.sequenceLength, self.nMonomers,
		self.sequenceMonomers, self.activeSequences from sequences.
		"""
		self.sequences = sequences
		self.nSequences, self.sequenceLength = sequences.shape
		self.nMonomers = np.amax(sequences) + 1 # [0 .. max(monomer #)]
		self.sequenceMonomers = np.empty(
			(self.nMonomers, self.nSequences, self.sequenceLength),
			dtype=bool
			)
		for monomerIndex in range(self.nMonomers):
			self.sequenceMonomers[monomerIndex, ...] = (sequences == monomerIndex)
		self.activeSequences = np.array(range(self.nSequences))

	def test_sum_monomers(self):
		sequences = np.array([
			[0, 0, 0, 0, P],
			[3, 3, 3, 3, P],
			[0, 1, 3, 0, 0]
			])
		self.makeSequenceMonomers(sequences)
		immediate = np.array([2, 0, 0, 1])

		indexes = np.full(self.activeSequences.shape, 0)
		tot = sum_monomers(self.sequenceMonomers[:, :], indexes, self.activeSequences)
		assert_equal(tot, immediate)

		tot2 = sum_monomers_reference_implementation(self.sequenceMonomers[:, :, 0],
			self.activeSequences)
		assert_equal(tot2, immediate)

	def test_partial_sum_monomers(self):
		sequences = np.array([
			[0, 0, 0, 0, P],
			[3, 3, 3, 3, P],
			[0, 0, 0, 1, 0],
			[2, 2, 2, 2, 2]
			])
		self.makeSequenceMonomers(sequences)
		self.activeSequences = np.array([0, 2]) # 2/4 sequences active
		currentStep = 2
		immediate = np.array([2, 0, 0, 0])

		indexes = np.full(self.activeSequences.shape, currentStep)
		tot = sum_monomers(
			self.sequenceMonomers[:, :],
			indexes,
			self.activeSequences)
		assert_equal(tot, immediate)

		tot2 = sum_monomers_reference_implementation(self.sequenceMonomers[:, :, currentStep],
			self.activeSequences)
		assert_equal(tot2, immediate)

	def test_polymerize_functionCall(self):
		sequences = np.array([
			[0, 1, 0, 1],
			[3, 1, 0, 2],
			[2, 2, 1, 1]
			])
		baseAmounts = np.array([9, 9, 9, 9])
		energy = 100

		rates = np.full(sequences.shape[0], 1)

		# Good calls test
		polymerize(sequences, baseAmounts, energy, np.random.RandomState(), rates)

	def test_polymerize_testEnergyLimited(self):
		sequences = np.array([
			[0, 1, 0, 1],
			[1, 1, 1, 1]
			])
		baseAmounts = np.array([9, 9])
		energy = 6

		rates = np.full(sequences.shape[0], 1)
		result = polymerize(sequences, baseAmounts, energy, np.random.RandomState(), rates)

		assert_equal(result.sequenceElongation, np.array([3, 3]))
		assert_equal(result.monomerUsages, np.array([2, 4]))
		self.assertEqual(6, result.nReactions)

	def test_polymerize_testFairness(self):
		sequences = np.array([
			[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
			[3] * 12,
			[2] * 12,
			[P] * 12])
		baseAmounts = np.array([11] * 4)
		baseAmountsOriginal = baseAmounts.copy()
		energy = 30

		rates = np.full(sequences.shape[0], 1)
		result = polymerize(sequences, baseAmounts, energy, np.random.RandomState(), rates)

		assert_equal(result.sequenceElongation, np.array([10, 9, 9, 0]))
		assert_equal(result.monomerUsages, np.array([3, 3, 11, 11]))
		self.assertEqual(28, result.nReactions)
		self.assertTrue((baseAmounts == baseAmountsOriginal).all(), 'should not modify monomerLimits')

	def test_polymerize_testAllAvailableBases(self):
		sequences = np.array([
			[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
			[3] * 12,
			[2] * 12,
			[P] * 12
			])
		baseAmounts = np.array([11] * 4)
		energy = 30

		rates = np.full(sequences.shape[0], 1)
		result = polymerize(sequences, baseAmounts, energy, np.random.RandomState(), rates)

		assert_equal(result.sequenceElongation, np.array([12,9,9,0]))
		assert_equal(result.monomerUsages, np.array([3,5,11,11]))
		self.assertEqual(30, result.nReactions)

	def test_polymerize_testVariableSequenceLength(self):
		sequences = np.array([
			[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
			[3] * 12,
			[2] * 12,
			[1, 1, 1] + [P] * 9
			])
		baseAmounts = np.array([30] * 4)
		energy = 50

		rates = np.full(sequences.shape[0], 1)
		result = polymerize(sequences, baseAmounts, energy, np.random.RandomState(), rates)

		assert_equal(result.sequenceElongation, np.array([12, 12, 12, 3]))

	def test_polymerize_testVariableElongation(self):
		sequences = np.array([
			[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
			[3] * 12,
			[2] * 12,
			[1, 1, 1] + [P] * 9])
		baseAmounts = np.array([30] * 4)
		energy = 20

		rates = np.array([4, 3, 2, 1])
		result = polymerize(
			sequences,
			baseAmounts,
			energy,
			np.random.RandomState(),
			rates,
			variable_elongation=True)

		assert_equal(result.sequenceElongation, np.array([8, 6, 4, 2]))

	def test_buildSequences(self):
		# Base case
		padding = np.empty((20, 10))
		padding.fill(P)
		allSequences =  np.hstack(
			(
				np.random.randint(3, size = (20, 10)), padding
				)
			).astype(np.int8, copy = False)
		sequenceIndexes = np.array([0, 5, 8])
		polymerizedLengths = np.array([0, 4, 9])
		elongationRate = 5
		rates = np.full(allSequences.shape[0], elongationRate)

		sequences = buildSequences(
			allSequences,
			sequenceIndexes,
			polymerizedLengths,
			rates)

		comparison_sequence = np.empty((3, elongationRate))
		comparison_sequence[0, :] = allSequences[0, 0:elongationRate]
		comparison_sequence[1, :] = allSequences[5, 4:(4+elongationRate)]
		comparison_sequence[2, :] = allSequences[8, 9:(9+elongationRate)]

		assert_equal(sequences, comparison_sequence)

		rates[-1] += 1000
		self.assertRaises(
			IndexError,
			buildSequences, allSequences, sequenceIndexes, polymerizedLengths, rates)

		# Un-padded case should throw exception
		allSequences = np.random.randint(3, size = (20, 10)).astype(np.int8, copy = False)
		sequenceIndexes = np.array([0, 5, 8])
		polymerizedLengths = np.array([0, 4, 9])
		elongationRate = 5

		self.assertRaises(
			TypeError,
			buildSequences, allSequences, sequenceIndexes, polymerizedLengths, elongationRate
			)

	def test_computeMassIncrease(self):
		sequences = np.random.randint(4, size = (20, 10)).astype(np.int8, copy = False)
		sequenceElongations = np.random.randint(10, size = (20,)).astype(np.int64, copy = False)
		monomerWeights = np.array([10., 4., 7., 1.], dtype = np.float64)

		massIncrease = computeMassIncrease(
			sequences,
			sequenceElongations,
			monomerWeights
			)

		comparison_massIncrease = np.empty(sequences.shape[0], dtype = np.float64)
		for i in range(sequences.shape[0]):
			comparison_massIncrease[i] = np.sum(monomerWeights[sequences[i][:sequenceElongations[i]]])

		assert_equal(massIncrease, comparison_massIncrease)
