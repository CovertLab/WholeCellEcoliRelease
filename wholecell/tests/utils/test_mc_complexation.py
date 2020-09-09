"""
Test polymerize_new.py
"""

from __future__ import absolute_import, division, print_function

from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

import numpy as np
import numpy.testing as npt

import unittest


class Test_mc_complexation(unittest.TestCase):

	def test_mccBuildMatrices(self):
		stoichMatrix = np.array([[-2, -1,  0,  0],
								 [ 1,  0, -1,  0],
								 [ 0,  1, -1,  0],
								 [ 0,  0,  1,  0],
								 [ 0,  0,  0, -1],
								 [ 0,  0,  0,  1]], dtype=np.int64)

		moleculeIndexes, overlappingReactions = mccBuildMatrices(stoichMatrix)

		npt.assert_equal(moleculeIndexes,
			np.array([[ 0,  1, -1],
					  [ 0,  2, -1],
					  [ 1,  2,  3],
					  [ 4,  5, -1]]))

		npt.assert_equal(overlappingReactions,
			np.array([[ 0,  1,  2],
					  [ 0,  1,  2],
					  [ 0,  1,  2],
					  [ 3, -1, -1]]))

	def test_mccFormComplexesWithPrebuiltMatrices(self):
		stoichMatrix = np.array([[-2, -1,  0,  0],
								 [ 1,  0, -1,  0],
								 [ 0,  1, -1,  0],
								 [ 0,  0,  1,  0],
								 [ 0,  0,  0, -1],
								 [ 0,  0,  0,  1]], dtype=np.int64)
		seed = 2
		test_moleculeCounts = np.array([[0, 6, 8, 0, 0, 0],
									[10, 0, 0, 0, 0, 0],
									[0, 0, 0, 3, 3, 0]], dtype=np.int64)

		expected_updatedMoleculeCounts = np.array([[0, 0, 2, 6, 0, 0],
									[0, 0, 1, 3, 0, 0],
									[0, 0, 0, 3, 0, 3]], dtype=np.int64)

		for i,test_moleculeCount in enumerate(test_moleculeCounts):
			prebuiltMatrices = mccBuildMatrices(stoichMatrix)

			updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
				test_moleculeCount,
				seed,
				stoichMatrix,
				*prebuiltMatrices
				)

			npt.assert_equal(updatedMoleculeCounts, expected_updatedMoleculeCounts[i])
