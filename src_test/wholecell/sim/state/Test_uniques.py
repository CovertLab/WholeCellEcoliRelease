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
import wholecell.sim.state.uniques

class Test_uniques(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.kb = type("", (), {})()
		self.kb.molecules = [{"id": "enz1", "mass": 1.0,"uniqueAttrs": None},
							{"id": "enz2", "mass": 2.0, "uniqueAttrs": None},
							{"id": "enz3", "mass": 3.0, "uniqueAttrs": ["attr1", "attr2", "attr3"]},
							{"id": "enz4", "mass": 4.0, "uniqueAttrs": ["attr4_1", "attr4_2"]},
							{"id": "enz5", "mass": 5.0, "uniqueAttrs": None}]
		self.kb.compartments = [{"id": "c"}, {"id": "e"}, {"id": "m"}]

		self.mc = wholecell.sim.state.uniques.MoleculesContainer()
		self.mc.initialize(self.kb)

	def tearDown(self):
		pass

	@noseAttrib.attr('uniqueTest')
	def test_countsBulk(self):
		mol = self.mc.molecule("enz1", "c")
		self.assertEqual(mol.countsBulk(), -1.0)

	@noseAttrib.attr('uniqueTest')
	def test_countsBulkIs(self):
		mol = self.mc.molecule("enz1", "c")
		mol.countsBulkIs(5)
		self.assertEqual(mol.countsBulk(), 5.0)
		mol.countsBulkIs(2)
		self.assertEqual(mol.countsBulk(), 2.0)

	@noseAttrib.attr('uniqueTest')
	def test_countsBulkInc(self):
		mol = self.mc.molecule("enz1", "c")
		mol.countsBulkIs(0)
		mol.countsBulkInc(1)
		self.assertEqual(mol.countsBulk(), 1.0)
		mol.countsBulkInc(2)
		self.assertEqual(mol.countsBulk(), 3.0)

	@noseAttrib.attr('uniqueTest')
	def test_countsBulkDec(self):
		mol = self.mc.molecule("enz1", "c")
		mol.countsBulkIs(3)
		mol.countsBulkDec(1)
		self.assertEqual(mol.countsBulk(), 2.0)
		mol.countsBulkDec(2)
		self.assertEqual(mol.countsBulk(), 0.0)
