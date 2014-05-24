#!/usr/bin/env python

"""
Test polymerize_new.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 1/22/2013
"""
from wholecell.utils.polymerize_new import polymerize

import numpy as np

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

class Test_polymerize(unittest.TestCase):
	
	@noseAttrib.attr('polymerizeNewTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_functionCall(self):
		sequences = np.array([[0,1,0,1],[3,1,0,2],[2,2,1,1]])
		baseAmounts = np.array([9,9,9,9])
		energy = 100
		
		# Good calls test
		polymerize(sequences, baseAmounts, energy)


	@noseAttrib.attr('polymerizeNewTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_testEnergyLimited(self):
		sequences = np.array([[0,1,0,1],[1,1,1,1]])
		baseAmounts = np.array([9,9])
		energy = 6
		
		np.random.seed(1)
		progress, baseCosts, energyCost = polymerize(sequences, baseAmounts, energy)
		np.random.seed()
		
		self.assertTrue((np.array([3, 3]) == progress).all())
		self.assertTrue((np.array([2, 4]) == baseCosts).all())
		self.assertEqual(6, energyCost)
	
	
	# @noseAttrib.attr('polymerizeNewTest')
	# @noseAttrib.attr('smalltest')
	# def test_polymerize_testFairness(self):
	# 	sequences = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
	# 							[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
	# 							[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
	# 							[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
	# 	baseAmounts = np.array([11,11,11,11])
	# 	bases = np.array([0,1,2,3])
	# 	basePadValue = ' '
	# 	energy = 30
	# 	energyCostPerBase = 1
		
	# 	np.random.seed(1)
	# 	progress, baseAmounts, baseCosts, energy, energyCost = polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
	# 	np.random.seed()
		
	# 	self.assertTrue((np.array([10,9,9,0]) == progress).all())
	# 	self.assertTrue((np.array([8,8,0,0]) == baseAmounts).all())
	# 	self.assertTrue((np.array([3,3,11,11]) == baseCosts).all())
	# 	self.assertEqual(2, energy)
	# 	self.assertEqual(28, energyCost)
	
	
	# @noseAttrib.attr('polymerizeNewTest')
	# @noseAttrib.attr('smalltest')
	# def test_polymerize_testAllAvailableBases(self):
	# 	sequences = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
	# 							[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
	# 							[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
	# 							[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
	# 	baseAmounts = np.array([11,11,11,11])
	# 	bases = np.array([0,1,2,3])
	# 	basePadValue = ' '
	# 	energy = 30
	# 	energyCostPerBase = 1
		
	# 	np.random.seed(1)
	# 	progress, baseAmounts, baseCosts, energy, energyCost = polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
	# 	np.random.seed()
		
	# 	self.assertTrue((np.array([12,9,9,0]) == progress).all())
	# 	self.assertTrue((np.array([8,6,0,0]) == baseAmounts).all())
	# 	self.assertTrue((np.array([3,5,11,11]) == baseCosts).all())
	# 	self.assertEqual(0, energy)
	# 	self.assertEqual(30, energyCost)
				
	
	
	# @noseAttrib.attr('polymerizeNewTest')
	# @noseAttrib.attr('smalltest')
	# def test_polymerize_testNoElongationOfFullyPolymerizedSequences(self):
		
		
	# 	# Polymerize with saturating tRNA, limiting energy, and 0th index
	# 	# sequence that has already been fully polymerized
	# 	elngSeqs = np.array([['0', '0'],
	# 							['24', '21'],
	# 							['29', '31']])
	# 	aminoacylatedTRNAs = 50 * np.ones((36,), dtype = np.int)
	# 	TRNAs = []
	# 	for i in range(1,36+1):
	# 		TRNAs.append(str(i))
	# 	TRNAs = np.array(TRNAs)
	# 	basePadValue = ' '
	# 	energy = 3
	# 	energyCostPerBase = 2
		
	# 	np.random.seed(3)
	# 	translationProgress, newAminoacylatedTRNAs, aminoacylatedTRNACosts, newEnergy, energyCost = polymerize(elngSeqs, aminoacylatedTRNAs, TRNAs, basePadValue, energy, energyCostPerBase)
	# 	np.random.seed()
		
		
	# 	# No elongation of fully polymerized sequences
	# 	self.assertEqual(0, translationProgress[np.where(elngSeqs[:, 0] == '0')[0][0]])
		
	# 	# Base costs correctly calculated
	# 	calculatedCounts = np.zeros((elngSeqs.shape[0],), dtype = np.int)
	# 	for i in range(elngSeqs.shape[1]):
	# 		pass
		
	# 	self.assertTrue((aminoacylatedTRNACosts == aminoacylatedTRNAs-newAminoacylatedTRNAs).all())
	# 	#assertEqual(aminoacylatedTRNACosts, tmp);
		
	# 	self.assertTrue((aminoacylatedTRNACosts == (aminoacylatedTRNAs - newAminoacylatedTRNAs)).all())
		
	# 	# Energy costs correctly calculated
	# 	self.assertTrue((energyCost == energy - newEnergy).all())
	# 	self.assertTrue((energyCost == energyCostPerBase * np.sum(aminoacylatedTRNACosts)).all())
	
	
	# @noseAttrib.attr('polymerizeNewTest')
	# @noseAttrib.attr('smalltest')
	# def test_polymerize_testVariableSequenceLength(self):
	# 	sequences = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 1],
	# 							[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
	# 							[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
	# 							[1, 1, 1, ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']])
	# 	baseAmounts = np.array([30,30,30,30])
	# 	bases = np.array([0,1,2,3])
	# 	basePadValue = ' '
	# 	energy = 50
	# 	energyCostPerBase = 1
		
	# 	np.random.seed(1)
	# 	progress, baseAmounts, baseCosts, energy, energyCost = polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
	# 	np.random.seed()
		
	# 	self.assertTrue((np.array([12,12,12,3]) == progress).all())
