#!/usr/bin/env python

"""
Test polymerize.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 1/22/2013
"""
import wholecell.utils.polymerize as p

import numpy as np
import pdb

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

class Test_polymerize(unittest.TestCase):
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_functionCall(self):
		sequences = np.matrix([['A','T','A','A'],['G','T','A','C'],['C','C','T','T']])
		baseAmounts = np.array([9,9,9,9])
		bases = np.array(['A','C','G','T'])
		basePadValue = ' '
		energy = 0
		energyCostPerBase = 0
		
		# Good calls test
		p.polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_validation_basic(self):
		sequences = np.matrix([['A','T','A','A'],['G','T','A','C'],['C','C','T','T']])
		baseAmounts = np.array([9,9,9,9])
		bases = np.array(['A','C','G','T'])
		basePadValue = ' '
		energy = 0
		energyCostPerBase = 0
		
		with self.assertRaises(p.polymerizeException) as context:
			sequences_bad = 'foobar'
			p.polymerize(sequences_bad, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'sequences must be a numpy matrix!\n')
		
		with self.assertRaises(p.polymerizeException) as context:
			baseAmounts_bad = 'foobar'
			p.polymerize(sequences, baseAmounts_bad, bases, basePadValue, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'baseAmounts must be an ndarray!\n')
		
		with self.assertRaises(p.polymerizeException) as context:
			bases_bad = 'foobar'
			p.polymerize(sequences, baseAmounts, bases_bad, basePadValue, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'bases must be an ndarray!\n')
		
		with self.assertRaises(p.polymerizeException) as context:
			basePadValue_bad = True
			p.polymerize(sequences, baseAmounts, bases, basePadValue_bad, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'basePadValue must be a string!\n')
	
			
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_validation_further(self):
		sequences = np.matrix([['A','T','A','A'],['G','T','A','C'],['C','C','T','T']])
		baseAmounts = np.array([9,9,9,9])
		bases = np.array(['A','C','G','T'])
		basePadValue = ' '
		energy = 0
		energyCostPerBase = 0
		
		with self.assertRaises(p.polymerizeException) as context:
			sequences_bad = np.matrix([], dtype = '|S1')
			p.polymerize(sequences_bad, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'sequences has no data!\n')
		
		with self.assertRaises(p.polymerizeException) as context:
			sequences_bad = np.matrix([['1','2','3','4'],['G','T','A','C'],['C','C','T','T']])
			p.polymerize(sequences_bad, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'sequences uses incorrect alphabet!\n')
		
		with self.assertRaises(p.polymerizeException) as context:
			basePadValue_bad = 'A'
			p.polymerize(sequences, baseAmounts, bases, basePadValue_bad, energy, energyCostPerBase)
		self.assertEqual(context.exception.message, 'basePadValue has the same value as a base!\n')
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testEnergyLimited(self):
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([9,9])
		bases = np.array(['A','B'])
		basePadValue = ' '
		energy = 6
		energyCostPerBase = 2
		
		np.random.seed(1)
		progress, baseAmounts, baseCosts, energy, energyCost = p.polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		np.random.seed()
		
		self.assertTrue((np.array([2,1]) == progress).all())
		self.assertTrue((np.array([8,7]) == baseAmounts).all())
		self.assertTrue((np.array([1,2]) == baseCosts).all())
		self.assertEqual(0, energy)
		self.assertEqual(6, energyCost)
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testFairness(self):
		sequences = np.matrix([['A', 'B', 'C', 'D', 'A', 'B', 'C', 'D', 'A', 'B', 'C', 'D'],
								['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D'],
								['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
								[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
		baseAmounts = np.array([11,11,11,11])
		bases = np.array(['A','B','C','D'])
		basePadValue = ' '
		energy = 30
		energyCostPerBase = 1
		
		np.random.seed(1)
		progress, baseAmounts, baseCosts, energy, energyCost = p.polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		np.random.seed()
		
		self.assertTrue((np.array([10,9,9,0]) == progress).all())
		self.assertTrue((np.array([8,8,0,0]) == baseAmounts).all())
		self.assertTrue((np.array([3,3,11,11]) == baseCosts).all())
		self.assertEqual(2, energy)
		self.assertEqual(28, energyCost)
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testAllAvailableBases(self):
		sequences = np.matrix([['A', 'B', 'C', 'D', 'A', 'B', 'C', 'D', 'A', 'B', 'B', 'B'],
								['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D'],
								['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
								[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
		baseAmounts = np.array([11,11,11,11])
		bases = np.array(['A','B','C','D'])
		basePadValue = ' '
		energy = 30
		energyCostPerBase = 1
		
		np.random.seed(1)
		progress, baseAmounts, baseCosts, energy, energyCost = p.polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		np.random.seed()
		
		self.assertTrue((np.array([12,9,9,0]) == progress).all())
		self.assertTrue((np.array([8,6,0,0]) == baseAmounts).all())
		self.assertTrue((np.array([3,5,11,11]) == baseCosts).all())
		self.assertEqual(0, energy)
		self.assertEqual(30, energyCost)
				
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testNoElongationOfFullyPolymerizedSequences(self):
		
		
		# Polymerize with saturating tRNA, limiting energy, and 0th index
		# sequence that has already been fully polymerized
		elngSeqs = np.matrix([['0', '0'],
								['24', '21'],
								['29', '31']])
		aminoacylatedTRNAs = 50 * np.ones((36,), dtype = np.int)
		TRNAs = []
		for i in range(1,36+1):
			TRNAs.append(str(i))
		TRNAs = np.array(TRNAs)
		basePadValue = ' '
		energy = 3
		energyCostPerBase = 2
		
		np.random.seed(3)
		translationProgress, newAminoacylatedTRNAs, aminoacylatedTRNACosts, newEnergy, energyCost = p.polymerize(elngSeqs, aminoacylatedTRNAs, TRNAs, basePadValue, energy, energyCostPerBase)
		np.random.seed()
		
		
		# No elongation of fully polymerized sequences
		self.assertEqual(0, translationProgress[np.where(elngSeqs[:, 0] == '0')[0][0]])
		
		# Base costs correctly calculated
		calculatedCounts = np.zeros((elngSeqs.shape[0],), dtype = np.int)
		for i in range(elngSeqs.shape[1]):
			pass
		
		self.assertTrue((aminoacylatedTRNACosts == aminoacylatedTRNAs-newAminoacylatedTRNAs).all())
		#assertEqual(aminoacylatedTRNACosts, tmp);
		
		self.assertTrue((aminoacylatedTRNACosts == (aminoacylatedTRNAs - newAminoacylatedTRNAs)).all())
		
		# Energy costs correctly calculated
		self.assertTrue((energyCost == energy - newEnergy).all())
		self.assertTrue((energyCost == energyCostPerBase * np.sum(aminoacylatedTRNACosts)).all())
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testVariableSequenceLength(self):
		sequences = np.matrix([['A', 'B', 'C', 'D', 'A', 'B', 'C', 'D', 'A', 'B', 'B', 'B'],
								['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D'],
								['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
								['B', 'B', 'B', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
		baseAmounts = np.array([30,30,30,30])
		bases = np.array(['A','B','C','D'])
		basePadValue = ' '
		energy = 50
		energyCostPerBase = 1
		
		np.random.seed(1)
		progress, baseAmounts, baseCosts, energy, energyCost = p.polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
		np.random.seed()
		
		self.assertTrue((np.array([12,12,12,3]) == progress).all())
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_calculateElongationLimits_stupidCaseProbablyUseless(self):
		sequences = np.matrix([['A', 'B', 'C', 'D', 'A', 'B', 'C', 'D', 'A', 'B', 'C', 'D'],
								['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D'],
								['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
								[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
		baseAmounts = np.array([11,11,11,11])
		bases = np.array(['A','B','C','D'])
		basePadValue = ' '
		energy = 30
		energyCostPerBase = 1		
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		self.assertEqual(elongation, 9)
		self.assertTrue((baseUsage == np.array([3,2,11,11])).all())
		self.assertTrue((limitingBases == np.array([2,3])).all())
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_calculateElongationLimits_enoughEnergyAndBases(self):
		# Enough energy and bases for multiple steps
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([9,9])
		bases = np.array(['A','B'])
		energy = 20
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 4
		test_baseUsage = np.array([2,6])
		test_limitingBases = np.array([])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_calculateElongationLimits_limitedEnergyAndBases(self):
		# Enough energy for one step and bases for multiple steps
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([9,9])
		bases = np.array(['A','B'])
		energy = 6
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 1
		test_baseUsage = np.array([1,1])
		test_limitingBases = np.array([])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
	
		
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_calculateElongationLimits_noEnergyEnoughBases(self):
		# Not enough energy for one step but enough bases for multiple steps
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([9,9])
		bases = np.array(['A','B'])
		energy = 2
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 0
		test_baseUsage = np.array([0,0])
		test_limitingBases = np.array([])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
		
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_calculateElongationLimits_enoughEnergyNoBases(self):
		# Enough energy for multiple steps but not enough bases for one step
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([0,0])
		bases = np.array(['A','B'])
		energy = 20
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 0
		test_baseUsage = np.array([0,0])
		test_limitingBases = np.array([0,1])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
		
	
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_calculateElongationLimits_enoughEnergyLimitedBases(self):
		# Enough energy for multiple steps but not enough bases for one step
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([1,0])
		bases = np.array(['A','B'])
		energy = 20
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 0
		test_baseUsage = np.array([0,0])
		test_limitingBases = np.array([1])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
		
		# Enough energy for multiple steps but not enough bases for one step
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([0,1])
		bases = np.array(['A','B'])
		energy = 20
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 0
		test_baseUsage = np.array([0,0])
		test_limitingBases = np.array([0])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
		
		# Enough energy for multiple steps but enough bases for a couple steps
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		baseAmounts = np.array([2,3])
		bases = np.array(['A','B'])
		energy = 20
		energyCostPerBase = 2
		
		elongation, baseUsage, limitingBases = p.calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
		
		test_elongation = 2
		test_baseUsage = np.array([1,3])
		test_limitingBases = np.array([1])
		
		self.assertEqual(test_elongation, elongation)
		self.assertTrue((test_baseUsage == baseUsage).all())
		self.assertTrue((test_limitingBases == limitingBases).all())
	
		
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_countBases(self):
		# Simple alphabet
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		bases = np.array(['A','B'])
		
		countBasesOutput = p.countBases(sequences, bases)
		
		test_countBasesOutput = np.matrix([[0,1,0,1,0], [0,1,2,1,2]])
		self.assertTrue((countBasesOutput == test_countBasesOutput).all())
		
		# Number alphabet
		sequences = np.matrix([['1','1','4','3','1'],['3','2','4','4','1'],['2','2','2','1','4']])
		bases = np.array(['1','2','3',4])
		
		countBasesOutput = p.countBases(sequences, bases)
		
		test_countBasesOutput = np.matrix([[0,1,1,0,1,2], [0,1,2,1,0,0], [0,1,0,0,1,0], [0,0,0,2,1,1]])
		self.assertTrue((countBasesOutput == test_countBasesOutput).all())
	
		
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_lengths(self):
		# No pad values
		sequences = np.matrix([['A','B','A','B'],['B','B','B','B']])
		padValue = ' '
		lengthsOutput = p.lengths(sequences, padValue)
		test_lengthsOutput = np.array([4,4])
		self.assertTrue((lengthsOutput == test_lengthsOutput).all())
		
		# With pad value
		sequences = np.matrix([['A','B','A',' '],['B','B','B','B']])
		padValue = ' '
		lengthsOutput = p.lengths(sequences, padValue)
		test_lengthsOutput = np.array([3,4])
		self.assertTrue((lengthsOutput == test_lengthsOutput).all())
		
		# With pad values
		sequences = np.matrix([['A','B','A',' '],['B','B','B',' '],['A','A','A','B']])
		padValue = ' '
		lengthsOutput = p.lengths(sequences, padValue)
		test_lengthsOutput = np.array([3,3,4])
		self.assertTrue((lengthsOutput == test_lengthsOutput).all())
	