#!/usr/bin/env python

"""
Test polymerize_new.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 6/23/2014
"""
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

import numpy as np

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

class Test_mc_complexation(unittest.TestCase):
	

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

	@noseAttrib.attr('complexationTest')
	@noseAttrib.attr('smalltest')
	def test_mccBuildMatrices(self):
		stoichMatrix = np.empty((2,1)).astype(np.int64, copy = False)
		stoichMatrix[0] = -2
		stoichMatrix[1] = 1
		stoichMatrix = np.array([[-2, -1,  0,  0],
								 [ 1,  0, -1,  0],
								 [ 0,  1, -1,  0],
								 [ 0,  0,  1,  0],
								 [ 0,  0,  0, -1],
								 [ 0,  0,  0,  1]], dtype=np.int64)

		moleculeIndexes, overlappingReactions = mccBuildMatrices(stoichMatrix)

		self.assertTrue(
			np.all(
				moleculeIndexes == np.array([[ 0,  1, -1],
											[ 0,  2, -1],
											[ 1,  2,  3],
											[ 4,  5, -1]])
				)
			)

		self.assertTrue(
			np.all(
				overlappingReactions == np.array([[ 0,  1],
											[ 0,  1],
											[ 0,  1],
											[ 3,  -1]])
				)
			)
