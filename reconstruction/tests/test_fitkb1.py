"""
Test fitkb1.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/23/2014
"""

from __future__ import division

import unittest
import warnings

import cPickle
import os

import wholecell.utils.constants
from reconstruction.ecoli.fitkb1 import totalCountFromMassesAndRatios, proteinDistributionFrommRNA, mRNADistributionFromProtein
import numpy as np
from wholecell.utils import units

import nose.plugins.attrib as noseAttrib

class Test_fitkb1(unittest.TestCase):

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

	@noseAttrib.attr('smalltest')
	@noseAttrib.attr('fitkb1test')
	def test_totalCountFromMassesAndRatios(self):
		# Test normal call
		totalMass = 10.
		individualMasses = np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.5, 0.25, 0.25])
		count = totalCountFromMassesAndRatios(totalMass, individualMasses, distribution)
		self.assertEqual(count, 80.)

		# Test normal call with units
		totalMass = 10. * units.fg
		individualMasses = units.fg * np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.5, 0.25, 0.25])
		count = totalCountFromMassesAndRatios(totalMass, individualMasses, distribution)
		count.checkNoUnit()
		self.assertEqual(count, 80.)

		# Test assertion in function
		totalMass = 10.
		individualMasses = np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.25, 0.25, 0.25])
		self.assertRaises(AssertionError, totalCountFromMassesAndRatios, totalMass, individualMasses, distribution)

	@noseAttrib.attr('smalltest')
	@noseAttrib.attr('fitkb1test')
	def test_proteinDistributionFrommRNA(self):
		# Test normal call
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = np.array([1, 2, 3])
		proteinDist = proteinDistributionFrommRNA(distribution_mRNA, netLossRate)
		distributionUnnormed = 1 / netLossRate * distribution_mRNA
		self.assertEqual(proteinDist.tolist(), (distributionUnnormed / units.sum(distributionUnnormed)).tolist())

		# Test normal call with units
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = proteinDistributionFrommRNA(distribution_mRNA, netLossRate)
		proteinDist.checkNoUnit()
		distributionUnnormed = 1 / netLossRate * distribution_mRNA
		self.assertEqual(proteinDist.asNumber().tolist(), (distributionUnnormed / units.sum(distributionUnnormed)).asNumber().tolist())

		# Test assertion in function
		distribution_mRNA = np.array([0.25, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		self.assertRaises(AssertionError, proteinDistributionFrommRNA, distribution_mRNA, netLossRate)

	@noseAttrib.attr('smalltest')
	@noseAttrib.attr('fitkb1test')
	def test_mRNADistributionFromProtein(self):
		# Test normal call
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = np.array([1, 2, 3])
		proteinDist = mRNADistributionFromProtein(distribution_mRNA, netLossRate)
		distributionUnnormed = netLossRate * distribution_mRNA
		self.assertEqual(proteinDist.tolist(), (distributionUnnormed / units.sum(distributionUnnormed)).tolist())

		# Test normal call with units
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = mRNADistributionFromProtein(distribution_mRNA, netLossRate)
		proteinDist.checkNoUnit()
		distributionUnnormed = netLossRate * distribution_mRNA
		self.assertEqual(proteinDist.asNumber().tolist(), (distributionUnnormed / units.sum(distributionUnnormed)).asNumber().tolist())

		# Test assertion in function
		distribution_mRNA = np.array([0.25, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		self.assertRaises(AssertionError, mRNADistributionFromProtein, distribution_mRNA, netLossRate)