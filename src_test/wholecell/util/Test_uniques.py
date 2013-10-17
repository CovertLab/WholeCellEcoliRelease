#!/usr/bin/env python

"""
Test uniques.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools

import numpy
import wholecell.util.randStream

class Test_uniques(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		cls.randStream = wholecell.util.randStream.randStream(seed = 1)

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.randStream.reset(1)

	def tearDown(self):
		pass


	@noseAttrib.attr('uniqueTest')
	def testProgram(self):
		r = self.randStream