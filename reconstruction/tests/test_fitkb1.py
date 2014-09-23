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
from reconstruction.ecoli.fitkb1 import totalCountFromMassesAndRatios
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
	def test_totalCountFromMassesAndRatios(self):
		totalMass = 10.
		individualMasses = np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.5, 0.25, 0.25])
		count = totalCountFromMassesAndRatios(totalMass, individualMasses, distribution)
		self.assertEqual(count, 80.)

		totalMass = 10. * units.fg
		individualMasses = units.fg * np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.5, 0.25, 0.25])
		count = totalCountFromMassesAndRatios(totalMass, individualMasses, distribution)
		count.checkNoUnit()
		self.assertEqual(count, 80.)

		totalMass = 10.
		individualMasses = np.array([0.1, 0.2, 0.1])
		distribution = np.array([0.25, 0.25, 0.25])
		self.assertRaises(AssertionError, totalCountFromMassesAndRatios, totalMass, individualMasses, distribution)
