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
