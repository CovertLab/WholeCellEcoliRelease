"""
Test fitkb1.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/23/2014
"""

from __future__ import absolute_import, division, print_function

import unittest

from reconstruction.ecoli.fit_sim_data_1 import (totalCountFromMassesAndRatios,
	proteinDistributionFrommRNA, mRNADistributionFromProtein,
	calculateMinPolymerizingEnzymeByProductDistribution,
	netLossRateFromDilutionAndDegradationProtein)

import numpy as np
from wholecell.utils import units


class Test_fitkb1(unittest.TestCase):

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
		self.assertFalse(units.hasUnit(count))
		self.assertEqual(count, 80.)

		# Test assertion in function
		totalMass = 10.
		individualMasses = np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.25, 0.25, 0.25])
		self.assertRaises(AssertionError, totalCountFromMassesAndRatios, totalMass, individualMasses, distribution)

	def test_proteinDistributionFrommRNA(self):
		# Test normal call
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = proteinDistributionFrommRNA(distribution_mRNA, np.ones(3) / 3, netLossRate)
		distributionUnnormed = 1 / netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		np.testing.assert_array_almost_equal_nulp(expectedDistribution.asNumber(), proteinDist)

		# Test normal call with units
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = proteinDistributionFrommRNA(distribution_mRNA, np.ones(3) / 3, netLossRate)
		distributionUnnormed = 1 / netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		np.testing.assert_array_almost_equal_nulp(expectedDistribution.asNumber(), proteinDist)

		# Test assertion in function
		distribution_mRNA = np.array([0.25, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		self.assertRaises(AssertionError, proteinDistributionFrommRNA, distribution_mRNA, np.ones(3) / 3, netLossRate)

	def test_mRNADistributionFromProtein(self):
		# Test normal call
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = mRNADistributionFromProtein(distribution_mRNA, np.ones(3) / 3, netLossRate)
		distributionUnnormed = netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		self.assertEqual(proteinDist.tolist(), expectedDistribution.asNumber().tolist())

		# Test normal call with units
		distribution_mRNA = np.array([0.5, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		proteinDist = mRNADistributionFromProtein(distribution_mRNA, np.ones(3) / 3, netLossRate)
		distributionUnnormed = netLossRate * distribution_mRNA
		expectedDistribution = (distributionUnnormed / units.sum(distributionUnnormed))
		expectedDistribution.normalize()
		expectedDistribution.checkNoUnit()
		self.assertEqual(proteinDist.tolist(), expectedDistribution.asNumber().tolist())

		# Test assertion in function
		distribution_mRNA = np.array([0.25, 0.25, 0.25])
		netLossRate = (1 / units.s) * np.array([1, 2, 3])
		self.assertRaises(AssertionError, mRNADistributionFromProtein, distribution_mRNA, np.ones(3) / 3, netLossRate)

	def test_calculateMinPolymerizingEnzymeByProductDistributionRNA(self):
		productLengths = units.aa * np.array([1000, 2000, 3000])
		elongationRates = units.aa / units.s * np.full(3, 50)
		netLossRate = (1 / units.s) * np.array([3, 2, 1])
		productCounts = np.array([5,1,10])

		nMin = calculateMinPolymerizingEnzymeByProductDistribution(productLengths, elongationRates, netLossRate, productCounts)
		nMin.checkNoUnit()
		self.assertEqual(nMin, 980)


	def test_netLossRateFromDilutionAndDegradationProtein(self):
		doublingTime = 60 * units.min
		degradationRates = (1 / units.s) * np.array([10, 20, 100])
		net = netLossRateFromDilutionAndDegradationProtein(doublingTime, degradationRates)
		net.asUnit(1/units.h) # Check units are 1/time
		self.assertEqual(
			net.asNumber(1/units.min).tolist(),
			((np.log(2) / doublingTime).asNumber(1/units.min) + degradationRates.asNumber(1/units.min)).tolist()
			)
