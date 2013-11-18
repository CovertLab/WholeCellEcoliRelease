#!/usr/bin/env python

"""
Test uniques.py

@author: Nick Ruggero, John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools

import numpy
# import wholecell.util.randStream
import wholecell.sim.state.uniques as wcUniques

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
							{"id": "enz5", "mass": 5.0, "uniqueAttrs": ["attr5_1"]}]
		self.kb.compartments = [{"id": "c"}, {"id": "e"}, {"id": "m"}]

		self.mc = wcUniques.MoleculesContainer()
		self.mc.initialize(self.kb)

	def tearDown(self):
		pass

	@noseAttrib.attr('uniqueTest')
	def test_moleculeNotUnique(self):
		try:
			mol = self.mc.molecule("enz1", "c")
		except:
			self.fail("Initalizing a molecule threw an error!")

	@noseAttrib.attr('uniqueTest')
	def test_countsBulk(self):
		mol = self.mc.molecule("enz1", "c")
		self.assertEqual(mol.countsBulk(), 0.0)

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

	@noseAttrib.attr('uniqueTest')
	def test_massSingle(self):
		mol = self.mc.molecule("enz2", "c")
		self.assertEqual(mol.massSingle(), 2.0)

	@noseAttrib.attr('uniqueTest')
	def test_massAllNoUnique(self):
		mol = self.mc.molecule("enz2", "c")
		mol.countsBulkIs(5)
		self.assertEqual(mol.massAll(), 10.0)

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		self.assertEqual(mol.countsUnique(), 1.0)
		self.assertEqual(newEnz3.attr1(), "A")
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")
		self.assertNotEqual(newEnz3.objects(), None)

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew_forNonUnique(self):
		mol = self.mc.molecule("enz1", "c")
		with self.assertRaises(wcUniques.uniqueException) as context:
			mol.uniqueNew({})
		self.assertEqual(context.exception.message, 'Attempting to create unique from object with no unique attributes!\n')

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew_missingCorrectAttr(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3 = mol.uniqueNew({"attr2" : "B", "attr3" : "C"})
		self.assertEqual(newEnz3.attr1(), None)
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")

		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr3" : "C"})
		self.assertEqual(newEnz3.attr1(), "A")
		self.assertEqual(newEnz3.attr2(), None)
		self.assertEqual(newEnz3.attr3(), "C")

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew_incorrectAttr(self):
		mol = self.mc.molecule("enz3", "c")
		with self.assertRaises(wcUniques.uniqueException) as context:
			newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C", "attr4" : "D"})
		self.assertEqual(context.exception.message, 'A specified attribute is not included in knoweldge base for this unique object!\n')

	@noseAttrib.attr('uniqueTest')
	def test_uniqueDel(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.uniqueDel(newEnz3)
		self.assertEqual(mol.countsUnique(), 0.0)

	@noseAttrib.attr('uniqueTest')
	def test_uniqueDel_objectNotCorrect(self):
		mol1 = self.mc.molecule("enz3", "c")
		mol2 = self.mc.molecule("enz4", "c")
		newEnz3 = mol1.uniqueNew()
		newEnz4 = mol2.uniqueNew()
		with self.assertRaises(wcUniques.uniqueException) as context:
			mol1.uniqueDel(newEnz4)
		self.assertEqual(context.exception.message, 'Unique object to delete does not match row in unique table!\n')

	@noseAttrib.attr('uniqueTest')
	def test_uniquesWithAttrs(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3_1 = mol.uniqueNew({"attr1" : "A", "attr2" : "A", "attr3" : "C"})
		newEnz3_2 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		newEnz3_3 = mol.uniqueNew({"attr1" : "B", "attr2" : "B", "attr3" : "C"})
		newEnz3_4 = mol.uniqueNew({"attr1" : "C", "attr2" : "B", "attr3" : "C"})

		L = mol.uniquesWithAttrs({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		self.assertEqual(len(L),1)
		self.assertEqual(id(L[0]), id(newEnz3_2))

		L = mol.uniquesWithAttrs({"attr1" : "A"})
		self.assertEqual(len(L),2)
		self.assertEqual(set([id(L[0]), id(L[1])]), set([id(newEnz3_1),id(newEnz3_2)]))

	@noseAttrib.attr('uniqueTest')
	def test_uniquesWithAttrs_incorrectAttr(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "A", "attr3" : "C"})
		with self.assertRaises(wcUniques.uniqueException) as context:
			mol.uniquesWithAttrs({"attr1" : "A", "attr4" : "B"})
		self.assertEqual(context.exception.message, 'A specified attribute is not included in knoweldge base for this unique object!\n')

	@noseAttrib.attr('uniqueTest')
	def test_makeGetter(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		self.assertEqual(newEnz3.attr1(), "A")
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")
		self.assertEqual(id(newEnz3), id(newEnz3.objects()))

	@noseAttrib.attr('uniqueTest')
	def test_makeSetter(self):
		mol = self.mc.molecule("enz3", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		newEnz3.attr1Is("Z")
		self.assertEqual(newEnz3.attr1(), "Z")
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")
		self.assertEqual(id(newEnz3), id(newEnz3.objects()))

	@noseAttrib.attr('uniqueTest')
	def test_dMassIs(self):
		mol = self.mc.molecule("enz3", "c")
		mol.countsBulkInc(1)
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.dMassIs(2)
		self.assertEqual(mol.massAll(), 8.0)

	@noseAttrib.attr('uniqueTest')
	def test_dMassInc(self):
		mol = self.mc.molecule("enz3", "c")
		mol.countsBulkInc(1)
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.dMassInc(2)
		self.assertEqual(mol.massAll(), 8.0)
		mol.dMassInc(1)
		self.assertEqual(mol.massAll(), 9.0)

	@noseAttrib.attr('uniqueTest')
	def test_dMassDec(self):
		mol = self.mc.molecule("enz3", "c")
		mol.countsBulkInc(1)
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.dMassIs(2)
		self.assertEqual(mol.massAll(), 8.0)
		mol.dMassDec(1)
		self.assertEqual(mol.massAll(), 7.0)

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_sameAsDefault(self):
		mol = self.mc.molecule("enz4", "c")
		e4 = mol.uniqueNew({"attr4_1" : "A", "attr4_2" : "B"})
		
		self.assertTrue('attr4_1' in dir(e4))
		self.assertTrue('attr4_1Is' in dir(e4))
		self.assertTrue('attr4_2' in dir(e4))
		self.assertTrue('attr4_2Is' in dir(e4))
		self.assertTrue('objects' in dir(e4))
		self.assertTrue('objectsIs' in dir(e4))

		self.assertEqual(e4.attr4_1(), "A")
		self.assertEqual(e4.attr4_2(), "B")
		self.assertEqual(id(e4.objects()), id(e4))
		e4.attr4_1Is("X")
		self.assertEqual(e4.attr4_1(), "X")
		e4.attr4_2Is("Y")
		self.assertEqual(e4.attr4_2(), "Y")

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_addNewFunction(self):
		mol = self.mc.molecule("enz5", "c")
		e5 = mol.uniqueNew({"attr5_1" : "A"})
		self.assertEqual(e5.test_function(), "Created new function")

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_keepsDefaultFunctions(self):
		mol = self.mc.molecule("enz5", "c")

		e5 = mol.uniqueNew({"attr5_1" : "A"})
		
		self.assertTrue('attr5_1' in dir(e5))
		self.assertTrue('attr5_1Is' in dir(e5))
		self.assertTrue('objects' in dir(e5))
		self.assertTrue('objectsIs' in dir(e5))

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_multipleUniqueMetaclasses(self):
		mol4 = self.mc.molecule("enz4", "c")
		e4 = mol4.uniqueNew({"attr4_1" : "A", "attr4_2" : "B"})
		mol5 = self.mc.molecule("enz5", "c")
		e5 = mol5.uniqueNew({"attr5_1" : "A"})

		self.assertTrue('attr4_1' in dir(e4))
		self.assertTrue('attr4_1Is' in dir(e4))
		self.assertTrue('attr4_2' in dir(e4))
		self.assertTrue('attr4_2Is' in dir(e4))
		self.assertTrue('objects' in dir(e4))
		self.assertTrue('objectsIs' in dir(e4))

		self.assertTrue('attr5_1' in dir(e5))
		self.assertTrue('attr5_1Is' in dir(e5))
		self.assertTrue('objects' in dir(e5))
		self.assertTrue('objectsIs' in dir(e5))
		self.assertTrue('test_function' in dir(e5))

class enz4Unique_same(object):
	registrationId = "enz4"
	__metaclass__ = wcUniques.MoleculeUniqueMeta

	def __init__(self, uniqueIdx):
		self._uniqueIdx = uniqueIdx

	def attr4_1(self):
		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_1"][self._uniqueIdx]

	def attr4_1Is(self, newVal):
		self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_1"][self._uniqueIdx] = newVal
		
	def attr4_2(self):
		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_2"][self._uniqueIdx]

	def attr4_2Is(self, newVal):
		self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_2"][self._uniqueIdx] = newVal

	def objects(self):
		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["objects"][self._uniqueIdx]

	def objectsIs(self, newVal):
		self._container._uniqueDict[self._molRowIdx][self._molColIdx]["objects"][self._uniqueIdx] = newVal

class enz5Unique_same(object):
	registrationId = "enz5"
	__metaclass__ = wcUniques.MoleculeUniqueMeta

	def test_function(self):
		return "Created new function"