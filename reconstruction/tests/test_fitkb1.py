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
from reconstruction.ecoli.fitkb1 import (totalCountFromMassesAndRatios,
proteinDistributionFrommRNA, mRNADistributionFromProtein,
calculateMinPolymerizingEnzymeByProductDistribution, netLossRateFromDilutionAndDegradation)

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
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = proteinDistributionFrommRNA(distribution_mRNA, netLossRate)
		distributionUnnormed = 1 / netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		self.assertEqual(proteinDist.tolist(), expectedDistribution.asNumber().tolist())

		# Test normal call with units
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = proteinDistributionFrommRNA(distribution_mRNA, netLossRate)
		distributionUnnormed = 1 / netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		self.assertEqual(proteinDist.tolist(), expectedDistribution.asNumber().tolist())

		# Test assertion in function
		distribution_mRNA = np.array([0.25, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		self.assertRaises(AssertionError, proteinDistributionFrommRNA, distribution_mRNA, netLossRate)

	@noseAttrib.attr('smalltest')
	@noseAttrib.attr('fitkb1test')
	def test_mRNADistributionFromProtein(self):
		# Test normal call
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = mRNADistributionFromProtein(distribution_mRNA, netLossRate)
		distributionUnnormed = netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		self.assertEqual(proteinDist.tolist(), expectedDistribution.asNumber().tolist())

		# Test normal call with units
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = mRNADistributionFromProtein(distribution_mRNA, netLossRate)
		distributionUnnormed = netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		self.assertEqual(proteinDist.tolist(), expectedDistribution.asNumber().tolist())

		# Test assertion in function
		distribution_mRNA = np.array([0.25, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		self.assertRaises(AssertionError, mRNADistributionFromProtein, distribution_mRNA, netLossRate)

	@noseAttrib.attr('smalltest')
	@noseAttrib.attr('fitkb1test')
	def test_calculateMinPolymerizingEnzymeByProductDistribution(self):
		productLengths = units.nt * np.array([1000, 2000, 3000])
		elongationRate = 50 * units.nt / units.s
		netLossRate = (1 / units.s) * np.array([3, 2, 1])
		productCounts = np.array([5,1,10])

		nMin = calculateMinPolymerizingEnzymeByProductDistribution(productLengths, elongationRate, netLossRate, productCounts)
		nMin.checkNoUnit()
		self.assertEqual(nMin, 980)


	@noseAttrib.attr('smalltest')
	@noseAttrib.attr('fitkb1test')
	def test_netLossRateFromDilutionAndDegradation(self):
		doublingTime = 60 * units.min
		degradationRates = (1 / units.s) * np.array([10, 20, 100])
		net = netLossRateFromDilutionAndDegradation(doublingTime, degradationRates)
		net.asUnit(1/units.h) # Check units are 1/time
		self.assertEqual(
			net.asNumber(1/units.min).tolist(),
			((np.log(2) / doublingTime).asNumber(1/units.min) + degradationRates.asNumber(1/units.min)).tolist()
			)
