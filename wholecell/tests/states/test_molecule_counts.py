#!/usr/bin/env python

"""
Test_MoleculeCounts.py

@author: Nick Ruggero, John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2013
"""

import unittest
import nose.plugins.attrib as noseAttrib
import cPickle
import os

import numpy as np
import wholecell.states.molecule_counts as wcMoleculeCounts

class Test_MoleculeCounts(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join("fixtures", "Simulation.cPickle"), "r"))
		self.kb = cPickle.load(open(os.path.join("fixtures","KnowledgeBase.cPickle"), "r"))

		self.bulkMolecules = self.sim.states['BulkMolecules']

		self.rand_counts = np.random.randint(0, 2000, size = self.bulkMolecules.countsBulk().shape)
		self.bulkMolecules.countsBulkIs(self.rand_counts)

	def tearDown(self):
		pass

	@noseAttrib.attr('smalltest')
	def test_CountsBulkViews(self):
		view1 = self.bulkMolecules.countsBulkViewNew()
		view2 = self.bulkMolecules.countsBulkViewNew([self.bulkMolecules._molIDs[0]
					 + '[' + self.bulkMolecules._compartments[0] + ']',])

		# Test accessing
		self.assertEqual(view1.countsBulk().tolist(),
			self.rand_counts.tolist())
		
		# Test modification within a view
		view1.countsBulkIs(view1.countsBulk() + 1)
		self.assertEqual((view1.countsBulk() - 1).tolist(),
			self.rand_counts.tolist())

		# Test modification across views
		self.assertEqual(view2.countsBulk() - 1, self.rand_counts[0,0])

	@noseAttrib.attr('smalltest')
	def test_CountsBulk(self):
		self.assertEqual(
			self.bulkMolecules.countsBulk().tolist(),
			self.rand_counts.tolist()
			)

	@noseAttrib.attr('smalltest')
	def test_PartitionIndexing(self):
		metabolite = 'ATP'
		compartment = 'c'

		partition_of_interest = self.bulkMolecules.partitions['Transcription']
		molIdx, cmpIdx = partition_of_interest._getIndex(
			'{}[{}]'.format(metabolite, compartment)
			)[1:]

		# Insure that _getIndex points to the right metabolite
		self.assertEqual(metabolite, partition_of_interest._molIDs[molIdx])
		self.assertEqual('merged', partition_of_interest._compartments[cmpIdx])

		# TODO: add more tests, more realistic tests cases


# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_moleculeNotUnique(self):
# 		try:
# 			mol = self.bulkMolecules.molecule("enz1[c]")
# 		except:
# 			self.fail("Initalizing a molecule threw an error!")

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_countBulk(self):
# 		mol = self.bulkMolecules.molecule("enz1[c]")
# 		self.assertEqual(mol.countBulk(), 0.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_countBulkIs(self):
# 		mol = self.bulkMolecules.molecule("enz1[c]")
# 		mol.countBulkIs(5)
# 		self.assertEqual(mol.countBulk(), 5.0)
# 		mol.countBulkIs(2)
# 		self.assertEqual(mol.countBulk(), 2.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_countBulkInc(self):
# 		mol = self.bulkMolecules.molecule("enz1[c]")
# 		mol.countBulkIs(0)
# 		mol.countBulkInc(1)
# 		self.assertEqual(mol.countBulk(), 1.0)
# 		mol.countBulkInc(2)
# 		self.assertEqual(mol.countBulk(), 3.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_countBulkDec(self):
# 		mol = self.bulkMolecules.molecule("enz1[c]")
# 		mol.countBulkIs(3)
# 		mol.countBulkDec(1)
# 		self.assertEqual(mol.countBulk(), 2.0)
# 		mol.countBulkDec(2)
# 		self.assertEqual(mol.countBulk(), 0.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_massSingle(self):
# 		mol = self.bulkMolecules.molecule("enz2[c]")
# 		self.assertEqual(mol.massSingle(), 2.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_massAllNoUnique(self):
# 		mol = self.bulkMolecules.molecule("enz2[c]")
# 		mol.countBulkIs(5)
# 		self.assertEqual(mol.massAll(), 10.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniqueNew(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		self.assertEqual(mol.countUnique(), 1.0)
# 		self.assertEqual(newEnz3.attr1(), "A")
# 		self.assertEqual(newEnz3.attr2(), "B")
# 		self.assertEqual(newEnz3.attr3(), "C")
# 		self.assertNotEqual(newEnz3.objects(), None)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniqueNew_forNonUnique(self):
# 		mol = self.bulkMolecules.molecule("enz1[c]")
# 		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
# 			mol.uniqueNew({})
# 		self.assertEqual(context.exception.message, 'Attempting to create unique from object with no unique attributes!\n')

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniqueNew_missingCorrectAttr(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3 = mol.uniqueNew({"attr2" : "B", "attr3" : "C"})
# 		self.assertEqual(newEnz3.attr1(), None)
# 		self.assertEqual(newEnz3.attr2(), "B")
# 		self.assertEqual(newEnz3.attr3(), "C")

# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr3" : "C"})
# 		self.assertEqual(newEnz3.attr1(), "A")
# 		self.assertEqual(newEnz3.attr2(), None)
# 		self.assertEqual(newEnz3.attr3(), "C")

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniqueNew_incorrectAttr(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
# 			newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C", "attr4" : "D"})
# 		self.assertEqual(context.exception.message, 'A specified attribute is not included in knoweldge base for this unique object!\n')

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniqueDel(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		mol.uniqueDel(newEnz3)
# 		self.assertEqual(mol.countUnique(), 0.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniqueDel_objectNotCorrect(self):
# 		mol1 = self.bulkMolecules.molecule("enz3[c]")
# 		mol2 = self.bulkMolecules.molecule("enz4[c]")
# 		newEnz3 = mol1.uniqueNew()
# 		newEnz4 = mol2.uniqueNew()
# 		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
# 			mol1.uniqueDel(newEnz4)
# 		self.assertEqual(context.exception.message, 'Unique object to delete does not match row in unique table!\n')

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniquesWithAttrs(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3_1 = mol.uniqueNew({"attr1" : "A", "attr2" : "A", "attr3" : "C"})
# 		newEnz3_2 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		newEnz3_3 = mol.uniqueNew({"attr1" : "B", "attr2" : "B", "attr3" : "C"})
# 		newEnz3_4 = mol.uniqueNew({"attr1" : "C", "attr2" : "B", "attr3" : "C"})

# 		L = mol.uniquesWithAttrs({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		self.assertEqual(len(L),1)
# 		self.assertEqual(id(L[0]), id(newEnz3_2))

# 		L = mol.uniquesWithAttrs({"attr1" : "A"})
# 		self.assertEqual(len(L),2)
# 		self.assertEqual(set([id(L[0]), id(L[1])]), set([id(newEnz3_1),id(newEnz3_2)]))

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_uniquesWithAttrs_incorrectAttr(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "A", "attr3" : "C"})
# 		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
# 			mol.uniquesWithAttrs({"attr1" : "A", "attr4" : "B"})
# 		self.assertEqual(context.exception.message, 'A specified attribute is not included in knoweldge base for this unique object!\n')

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_makeGetter(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		self.assertEqual(newEnz3.attr1(), "A")
# 		self.assertEqual(newEnz3.attr2(), "B")
# 		self.assertEqual(newEnz3.attr3(), "C")
# 		self.assertEqual(id(newEnz3), id(newEnz3.objects()))

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_makeSetter(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		newEnz3.attr1Is("Z")
# 		self.assertEqual(newEnz3.attr1(), "Z")
# 		self.assertEqual(newEnz3.attr2(), "B")
# 		self.assertEqual(newEnz3.attr3(), "C")
# 		self.assertEqual(id(newEnz3), id(newEnz3.objects()))

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_dMassIs(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		mol.countBulkInc(1)
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		mol.dMassIs(2)
# 		self.assertEqual(mol.massAll(), 8.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_dMassInc(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		mol.countBulkInc(1)
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		mol.dMassInc(2)
# 		self.assertEqual(mol.massAll(), 8.0)
# 		mol.dMassInc(1)
# 		self.assertEqual(mol.massAll(), 9.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_dMassDec(self):
# 		mol = self.bulkMolecules.molecule("enz3[c]")
# 		mol.countBulkInc(1)
# 		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
# 		mol.dMassIs(2)
# 		self.assertEqual(mol.massAll(), 8.0)
# 		mol.dMassDec(1)
# 		self.assertEqual(mol.massAll(), 7.0)

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_metaClass_sameAsDefault(self):
# 		mol = self.bulkMolecules.molecule("enz4[c]")
# 		e4 = mol.uniqueNew({"attr4_1" : "A", "attr4_2" : "B"})

# 		self.assertTrue(hasattr(e4, 'attr4_1'))
# 		self.assertTrue(hasattr(e4, 'attr4_1Is'))
# 		self.assertTrue(hasattr(e4, 'attr4_2'))
# 		self.assertTrue(hasattr(e4, 'attr4_2Is'))
# 		self.assertTrue(hasattr(e4, 'objects'))
# 		self.assertTrue(hasattr(e4, 'objectsIs'))

# 		self.assertEqual(e4.attr4_1(), "A")
# 		self.assertEqual(e4.attr4_2(), "B")
# 		self.assertEqual(id(e4.objects()), id(e4))
# 		e4.attr4_1Is("X")
# 		self.assertEqual(e4.attr4_1(), "X")
# 		e4.attr4_2Is("Y")
# 		self.assertEqual(e4.attr4_2(), "Y")

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_metaClass_addNewFunction(self):
# 		mol = self.bulkMolecules.molecule("enz5[c]")
# 		e5 = mol.uniqueNew({"attr5_1" : "A"})
# 		self.assertEqual(e5.test_function(), "Created new function")

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_metaClass_keepsDefaultFunctions(self):
# 		mol = self.bulkMolecules.molecule("enz5[c]")

# 		e5 = mol.uniqueNew({"attr5_1" : "A"})

# 		self.assertTrue(hasattr(e5, 'attr5_1'))
# 		self.assertTrue(hasattr(e5, 'attr5_1Is'))
# 		self.assertTrue(hasattr(e5, 'objects'))
# 		self.assertTrue(hasattr(e5, 'objectsIs'))

# # 	@noseAttrib.attr('brokensmalltest') #@noseAttrib.attr('smalltest')
# 	def test_metaClass_multipleUniqueMetaclasses(self):
# 		mol4 = self.bulkMolecules.molecule("enz4[c]")
# 		e4 = mol4.uniqueNew({"attr4_1" : "A", "attr4_2" : "B"})
# 		mol5 = self.bulkMolecules.molecule("enz5[c]")
# 		e5 = mol5.uniqueNew({"attr5_1" : "A"})

# 		self.assertTrue(hasattr(e4, 'attr4_1'))
# 		self.assertTrue(hasattr(e4, 'attr4_1Is'))
# 		self.assertTrue(hasattr(e4, 'attr4_2'))
# 		self.assertTrue(hasattr(e4, 'attr4_2Is'))
# 		self.assertTrue(hasattr(e4, 'objects'))
# 		self.assertTrue(hasattr(e4, 'objectsIs'))

# 		self.assertTrue(hasattr(e5, 'attr5_1'))
# 		self.assertTrue(hasattr(e5, 'attr5_1Is'))
# 		self.assertTrue(hasattr(e5, 'objects'))
# 		self.assertTrue(hasattr(e5, 'objectsIs'))
# 		self.assertTrue(hasattr(e5, 'test_function'))


# class enz4Unique_same(object):
# 	registrationId = "enz4"
# 	__metaclass__ = wcMoleculeCounts.MoleculeUniqueMeta

# 	def __init__(self, uniqueIdx):
# 		self._uniqueIdx = uniqueIdx

# 	def attr4_1(self):
# 		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_1"][self._uniqueIdx]

# 	def attr4_1Is(self, newVal):
# 		self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_1"][self._uniqueIdx] = newVal
		
# 	def attr4_2(self):
# 		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_2"][self._uniqueIdx]

# 	def attr4_2Is(self, newVal):
# 		self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_2"][self._uniqueIdx] = newVal

# 	def objects(self):
# 		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["objects"][self._uniqueIdx]

# 	def objectsIs(self, newVal):
# 		self._container._uniqueDict[self._molRowIdx][self._molColIdx]["objects"][self._uniqueIdx] = newVal


# class enz5Unique_same(object):
# 	registrationId = "enz5"
# 	__metaclass__ = wcMoleculeCounts.MoleculeUniqueMeta

# 	def test_function(self):
# 		return "Created new function"

