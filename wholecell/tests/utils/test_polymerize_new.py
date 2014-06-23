#!/usr/bin/env python

"""
Test polymerize_new.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 1/22/2013
"""
from wholecell.utils.polymerize_new import  buildSequences, polymerize, computeMassIncrease, PAD_VALUE

import numpy as np

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
		
		self.assertTrue((np.array([3, 3]) == progress).all())
		self.assertTrue((np.array([2, 4]) == baseCosts).all())
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
		
		self.assertTrue((np.array([10,9,9,0]) == progress).all())
		self.assertTrue((np.array([3,3,11,11]) == baseCosts).all())
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
		
		self.assertTrue((np.array([12,9,9,0]) == progress).all())
		self.assertTrue((np.array([3,5,11,11]) == baseCosts).all())
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
		
		self.assertTrue((np.array([12,12,12,3]) == progress).all())


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

		self.assertTrue(np.all(
			sequences == comparison_sequence
			))

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

		self.assertTrue(np.all(massIncrease == comparison_massIncrease))