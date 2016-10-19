"""
Test polymerize.py

	cd wcEcoli
	nosetests wholecell/tests/utils/test_polymerize.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 1/22/2013
"""

from wholecell.utils.polymerize import (buildSequences, polymerize,
	computeMassIncrease, PAD_VALUE, sum_monomers,
	sum_monomers_reference_implementation)

import numpy as np
from numpy.testing import assert_equal

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest


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
		self.sequenceMonomers = np.empty((self.nMonomers, self.nSequences, self.sequenceLength), dtype=np.bool)
		for monomerIndex in xrange(self.nMonomers):
			self.sequenceMonomers[monomerIndex, ...] = (sequences == monomerIndex)
		self.activeSequences = np.array(xrange(self.nSequences))

	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_sum_monomers(self):
		sequences = np.array([[0,0,0,0,PAD_VALUE], [3,3,3,3,PAD_VALUE], [0,1,3,0,0]])
		self.makeSequenceMonomers(sequences)
		expected = np.array([[2,3,4,6,7], [0,1,1,1,1], [0,0,0,0,0], [1,2,4,5,5]])
		tot = sum_monomers(self.sequenceMonomers, self.activeSequences, 0)
		assert_equal(tot, expected)

		tot2 = sum_monomers_reference_implementation(self.sequenceMonomers,
			self.activeSequences, 0)
		assert_equal(tot2, expected)

	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_partial_sum_monomers(self):
		sequences = np.array([[0,0, 0,0,PAD_VALUE], [3,3, 3,3,PAD_VALUE],
			[0,0, 0,1,0], [2,2, 2,2,2]])
		self.makeSequenceMonomers(sequences)
		self.activeSequences = np.array([0,2]) # 2/4 sequences active
		currentStep = 2
		expected = np.array([[2,3,4], [0,1,1], [0,0,0], [0,0,0]])
		tot = sum_monomers(self.sequenceMonomers, self.activeSequences, currentStep)
		assert_equal(tot, expected)

		tot2 = sum_monomers_reference_implementation(self.sequenceMonomers,
			self.activeSequences, currentStep)
		assert_equal(tot2, expected)

	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_polymerize_functionCall(self):
		sequences = np.array([[0,1,0,1],[3,1,0,2],[2,2,1,1]])
		baseAmounts = np.array([9,9,9,9])
		energy = 100
		
		# Good calls test
		polymerize(sequences, baseAmounts, energy, np.random.RandomState(0))


	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testEnergyLimited(self):
		sequences = np.array([[0,1,0,1],[1,1,1,1]])
		baseAmounts = np.array([9,9])
		energy = 6
		
		np.random.seed(1)
		progress, baseCosts, energyCost = polymerize(sequences, baseAmounts, energy, np.random.RandomState(0))
		np.random.seed()
		
		assert_equal(progress, np.array([3, 3]))
		assert_equal(baseCosts, np.array([2, 4]))
		self.assertEqual(6, energyCost)
	
	
	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testFairness(self):
		sequences = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
								[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
								[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
								[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]])
		baseAmounts = np.array([11,11,11,11])
		energy = 30
		
		np.random.seed(1)
		progress, baseCosts, energyCost = polymerize(sequences, baseAmounts, energy, np.random.RandomState(0))
		np.random.seed()
		
		assert_equal(progress, np.array([10,9,9,0]))
		assert_equal(baseCosts, np.array([3,3,11,11]))
		self.assertEqual(28, energyCost)
	
	
	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testAllAvailableBases(self):
		sequences = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
								[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
								[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
								[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]])
		baseAmounts = np.array([11,11,11,11])
		energy = 30
		
		np.random.seed(1)
		progress, baseCosts, energyCost = polymerize(sequences, baseAmounts, energy, np.random.RandomState(0))
		np.random.seed()
		
		assert_equal(progress, np.array([12,9,9,0]))
		assert_equal(baseCosts, np.array([3,5,11,11]))
		self.assertEqual(30, energyCost)
				
	
	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testVariableSequenceLength(self):
		sequences = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
								[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
								[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
								[1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1]])
		baseAmounts = np.array([30,30,30,30])
		energy = 50
		
		np.random.seed(1)
		progress, baseCosts, energyCost = polymerize(sequences, baseAmounts, energy, np.random.RandomState(0))
		np.random.seed()
		
		assert_equal(progress, np.array([12,12,12,3]))


	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_buildSequences(self):
		# Base case
		padding = np.empty((20,10))
		padding.fill(PAD_VALUE)
		allSequences =  np.hstack(
			(
				np.random.randint(3, size=(20,10)), padding
				)
			).astype(np.int8, copy=False)
		sequenceIndexes = np.array([0,5,8])
		polymerizedLengths = np.array([0,4,9])
		elngRate = 5
		
		sequences = buildSequences(
			allSequences,
			sequenceIndexes,
			polymerizedLengths,
			elngRate
			)

		comparison_sequence = np.empty((3,elngRate))
		comparison_sequence[0,:] = allSequences[0,0:elngRate]
		comparison_sequence[1,:] = allSequences[5,4:4+elngRate]
		comparison_sequence[2,:] = allSequences[8,9:9+elngRate]

		assert_equal(sequences, comparison_sequence)

		# Un-padded case should throw exception
		allSequences = np.random.randint(3, size=(20,10)).astype(np.int8, copy=False)
		sequenceIndexes = np.array([0,5,8])
		polymerizedLengths = np.array([0,4,9])
		elngRate = 5
		
		# TODO: Define an exeption subclass if you want to check the real message
		self.assertRaises(
			Exception,
			buildSequences, allSequences, sequenceIndexes, polymerizedLengths, elngRate
			)

	@noseAttrib.attr('polymerizeNew')
	@noseAttrib.attr('smalltest')
	def test_computeMassIncrease(self):
		sequences = np.random.randint(4, size=(20,10)).astype(np.int8, copy=False)
		sequenceElongations = np.random.randint(10, size = (20,)).astype(np.int64, copy=False)
		monomerWeights = np.array([10., 4., 7., 1.], dtype=np.float64)
		
		massIncrease = computeMassIncrease(
			sequences,
			sequenceElongations,
			monomerWeights
			)

		comparison_massIncrease = np.empty(sequences.shape[0], dtype=np.float64)
		for i in range(sequences.shape[0]):
			comparison_massIncrease[i] = np.sum(monomerWeights[sequences[i][:sequenceElongations[i]]])

		assert_equal(massIncrease, comparison_massIncrease)
