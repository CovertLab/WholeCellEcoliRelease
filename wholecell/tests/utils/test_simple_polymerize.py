#!/usr/bin/env python

"""
Test simple_polymerize.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/22/2014
"""

from __future__ import division

import unittest
import warnings

import numpy as np
from wholecell.utils.rand_stream import RandStream
from wholecell.utils.simple_polymerize import simplePolymerize

import nose.plugins.attrib as noseAttrib

class Test_simplePolymerize(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.randStream = RandStream(seed = 1)

	def tearDown(self):
		pass

	@noseAttrib.attr('mediumtest', 'process', 'polymerize')
	def test_flat_distribution(self):
		# Test behavior for identical monomers and synthesis probabilities
		nMonomerTypes = 4
		nMonomersEachType = 1000

		nTemplates = 20
		nEachMonomerOnTemplate = 4

		templateMonomerCounts = nEachMonomerOnTemplate * np.ones(
			(nTemplates, nMonomerTypes), np.int)

		synthesisProbabilities = 1/nTemplates * np.ones(nTemplates)

		self.assertTrue(np.abs(synthesisProbabilities.sum() - 1) < 1e-5)

		totalPolymerCounts = np.zeros(nTemplates, np.int)

		nRounds = 1000

		for i in xrange(nRounds):
			enzymaticLimitation = 1e9 # no limitation
			monomerCounts = nMonomersEachType * np.ones(nMonomerTypes, np.int)

			polymerCounts = simplePolymerize(
				templateMonomerCounts,
				enzymaticLimitation,
				monomerCounts,
				synthesisProbabilities,
				self.randStream
				)[1]

			totalPolymerCounts += polymerCounts
		
		observedSynthProb = totalPolymerCounts / totalPolymerCounts.sum()

		self.assertTrue(
			(np.abs(observedSynthProb - synthesisProbabilities) < 1e-2).all()
			)


	@noseAttrib.attr('mediumtest', 'process', 'polymerize')
	def test_nonflat_distribution(self):
		# Test behavior for varied synthesis probabilities
		nMonomerTypes = 4
		nMonomersEachType = 1000

		nTemplates = 20
		nEachMonomerOnTemplate = 4

		templateMonomerCounts = nEachMonomerOnTemplate * np.ones(
			(nTemplates, nMonomerTypes), np.int)

		# Ramp up
		synthesisProbabilities = 1.5**(np.arange(nTemplates) + 1.)
		synthesisProbabilities /= synthesisProbabilities.sum()

		self.assertTrue(np.abs(synthesisProbabilities.sum() - 1) < 1e-5)

		totalPolymerCounts = np.zeros(nTemplates, np.int)

		nRounds = 1000

		for i in xrange(nRounds):
			enzymaticLimitation = 1e9 # no limitation
			monomerCounts = nMonomersEachType * np.ones(nMonomerTypes, np.int)

			polymerCounts = simplePolymerize(
				templateMonomerCounts,
				enzymaticLimitation,
				monomerCounts,
				synthesisProbabilities,
				self.randStream
				)[1]

			totalPolymerCounts += polymerCounts
		
		observedSynthProb = totalPolymerCounts / totalPolymerCounts.sum()

		self.assertTrue(
			(np.abs(observedSynthProb - synthesisProbabilities) < 1e-2).all()
			)






		