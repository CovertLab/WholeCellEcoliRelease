#!/usr/bin/env python

"""
Test MoleculeCounts.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/5/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools

import numpy
import cPickle
import os
# import wholecell.util.randStream
import wholecell.sim.state.MoleculeCounts


class Test_MoleculeCounts(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = None
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		self.moleculeCounts = self.sim.getState('MoleculeCounts')

	def tearDown(self):
		pass

	@noseAttrib.attr('partitionTest')
	def test_Basic(self):
		import ipdb
		ipdb.set_trace()
		
		self.moleculeCounts.prepartition()
		self.moleculeCounts.partition()
		ipdb.set_trace()