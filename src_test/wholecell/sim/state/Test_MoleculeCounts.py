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

import numpy
import wholecell.sim.state.MoleculeCounts as wcMoleculeCounts

class Test_MoleculeCounts(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp_new(self):
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		self.mc = self.sim.getState('MoleculeCounts')

	def setUp(self):
		self.kb = type("KnowledgeBase", (object,), {'molecules':None})()
		self.kb.molecules = [
			{"id": "enz1:mature", "mass": 1.0,"uniqueAttrs": None},
			{"id": "enz2:mature", "mass": 2.0, "uniqueAttrs": None},
			{"id": "enz3:mature", "mass": 3.0, "uniqueAttrs": ["attr1", "attr2", "attr3"]},
			{"id": "enz4:mature", "mass": 4.0, "uniqueAttrs": ["attr4_1", "attr4_2"]},
			{"id": "enz5:mature", "mass": 5.0, "uniqueAttrs": ["attr5_1"]},
			{"id": "metA:mature", "mass": 1.0, "uniqueAttrs": None},
			{"id": "metB:mature", "mass": 1.0, "uniqueAttrs": None},
			{"id": "metC:mature", "mass": 1.0, "uniqueAttrs": None},
			{"id": "metD:mature", "mass": 1.0, "uniqueAttrs": None},
			{"id": "metE:mature", "mass": 1.0, "uniqueAttrs": None},
			{"id": "metF:mature", "mass": 1.0, "uniqueAttrs": None},
			{"id": "metG:mature", "mass": 1.0, "uniqueAttrs": None},
			]
		self.kb.compartments = [{"id": "c"}, {"id": "e"}, {"id": "m"}]

		self.enzIds = [mol['id'] for mol in self.kb.molecules if 'enz' in mol['id']]
		self.metIds = [mol['id'] for mol in self.kb.molecules if 'met' in mol['id']]

		self.metCounts = [10.,2.,5.,7.,20.,3.,7.]

		self.mc = wcMoleculeCounts.MoleculeCounts()
		self.mc.initialize_old(self.kb)

		# Create generic process for partition
		self.genericProcess = type("Process", (object,), {'meta':None})()
		self.genericProcess.meta = {"id": "genericProcess_id", "name": "genericProcess_name"}

		# Create some partitions, currently without request functions
		self.partition1 = self.mc.addPartition(self.genericProcess, [
			metId + '[c]' for metId in self.metIds], lambda partition: None)

		self.partition2 = self.mc.addPartition(self.genericProcess, [
			metId + '[c]' for metId in self.metIds], lambda partition: None)

		self.partition3 = self.mc.addPartition(self.genericProcess, [
			metId + '[c]' for metId in self.metIds], lambda partition: None)

		self.mc.allocate()

		self.mc.countsBulkIs(self.metCounts, self.metIds)

	def tearDown(self):
		pass

	@noseAttrib.attr('uniqueTest')
	def test_CountsBulkViews(self):
		view1 = self.mc.countsBulkViewNew(self.metIds)
		view2 = self.mc.countsBulkViewNew([self.metIds[0],])

		# Test accessing
		self.assertEqual(view1.countsBulk().tolist(),
			self.metCounts)
		
		view1.countsBulkIs(view1.countsBulk() + 1)

		# Test modification within a view
		self.assertEqual((view1.countsBulk() - 1).tolist(),
			self.metCounts)

		# Test modification across views
		self.assertEqual(view2.countsBulk() - 1, self.metCounts[0])


	@noseAttrib.attr('uniqueTest')
	def test_CountsBulkViewsPartition(self):
		view = self.partition1.countsBulkViewNew([self.metIds[0], ])

		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.mc.prepartition()
		self.mc.partition()

		self.assertEqual(view.countsBulk(), 1.)


	@noseAttrib.attr('uniqueTest')
	def test_CountsBulk(self):
		self.assertEqual(
			self.mc.countsBulk(self.metIds).tolist(),
			self.metCounts
			)


	@noseAttrib.attr('uniqueTest')
	def test_PartitionIndexing(self):
		metaboliteID = 'metA:mature'

		# Insure that _getIndex points to the right metabolite
		self.assertEqual(metaboliteID, self.partition1._wids[
			self.partition1._getIndex(metaboliteID)[0]
			])

		# TODO: add more tests, more realistic tests cases

	@noseAttrib.attr('uniqueTest')
	def test_relativeAllocation(self):
		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.mc.prepartition()
		self.mc.partition()

		self.assertEqual(self.partition1.countsBulk().flatten().tolist(), [1., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.partition2.countsBulk().flatten().tolist(), [1., 1., 2., 4., 0., 0., 0.])
		self.assertEqual(self.partition3.countsBulk().flatten().tolist(), [7., 0., 3., 2., 20., 0., 7.])

	@noseAttrib.attr('uniqueTest')
	def test_absoluteAllocation(self):
		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.partition2.isReqAbs = True # this is a hack

		self.mc.prepartition()
		self.mc.partition()
		
		self.assertEqual(self.partition1.countsBulk().flatten().tolist(), [0., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.partition2.countsBulk().flatten().tolist(), [5., 2., 2., 2., 0., 0., 0.])
		self.assertEqual(self.partition3.countsBulk().flatten().tolist(), [4., 0., 3., 5., 20., 0., 7.])

	@noseAttrib.attr('uniqueTest')
	def test_absoluteAllocation_withConflict(self):
		self.partition1.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.partition2.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.partition3.reqFunc = lambda partition: partition.countsBulkIs(numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.partition2.isReqAbs = True # this is a hack
		self.partition3.isReqAbs = True # this is a hack

		self.mc.prepartition()
		self.mc.partition()
		
		self.assertEqual(self.partition1.countsBulk().flatten().tolist(), [0., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.partition2.countsBulk().flatten().tolist(), [2., 1., 2., 2., 0., 0., 0.])
		self.assertEqual(self.partition3.countsBulk().flatten().tolist(), [8., 0., 3., 1., 2., 0., 2.])

	@noseAttrib.attr('uniqueTest')
	def test_moleculeNotUnique(self):
		try:
			mol = self.mc.molecule("enz1:mature", "c")
		except:
			self.fail("Initalizing a molecule threw an error!")

	@noseAttrib.attr('uniqueTest')
	def test_countBulk(self):
		mol = self.mc.molecule("enz1:mature", "c")
		self.assertEqual(mol.countBulk(), 0.0)

	@noseAttrib.attr('uniqueTest')
	def test_countBulkIs(self):
		mol = self.mc.molecule("enz1:mature", "c")
		mol.countBulkIs(5)
		self.assertEqual(mol.countBulk(), 5.0)
		mol.countBulkIs(2)
		self.assertEqual(mol.countBulk(), 2.0)

	@noseAttrib.attr('uniqueTest')
	def test_countBulkInc(self):
		mol = self.mc.molecule("enz1:mature", "c")
		mol.countBulkIs(0)
		mol.countBulkInc(1)
		self.assertEqual(mol.countBulk(), 1.0)
		mol.countBulkInc(2)
		self.assertEqual(mol.countBulk(), 3.0)

	@noseAttrib.attr('uniqueTest')
	def test_countBulkDec(self):
		mol = self.mc.molecule("enz1:mature", "c")
		mol.countBulkIs(3)
		mol.countBulkDec(1)
		self.assertEqual(mol.countBulk(), 2.0)
		mol.countBulkDec(2)
		self.assertEqual(mol.countBulk(), 0.0)

	@noseAttrib.attr('uniqueTest')
	def test_massSingle(self):
		mol = self.mc.molecule("enz2:mature", "c")
		self.assertEqual(mol.massSingle(), 2.0)

	@noseAttrib.attr('uniqueTest')
	def test_massAllNoUnique(self):
		mol = self.mc.molecule("enz2:mature", "c")
		mol.countBulkIs(5)
		self.assertEqual(mol.massAll(), 10.0)

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew(self):
		mol = self.mc.molecule("enz3:mature", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		self.assertEqual(mol.countUnique(), 1.0)
		self.assertEqual(newEnz3.attr1(), "A")
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")
		self.assertNotEqual(newEnz3.objects(), None)

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew_forNonUnique(self):
		mol = self.mc.molecule("enz1:mature", "c")
		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
			mol.uniqueNew({})
		self.assertEqual(context.exception.message, 'Attempting to create unique from object with no unique attributes!\n')

	@noseAttrib.attr('uniqueTest')
	def test_uniqueNew_missingCorrectAttr(self):
		mol = self.mc.molecule("enz3:mature", "c")
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
		mol = self.mc.molecule("enz3:mature", "c")
		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
			newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C", "attr4" : "D"})
		self.assertEqual(context.exception.message, 'A specified attribute is not included in knoweldge base for this unique object!\n')

	@noseAttrib.attr('uniqueTest')
	def test_uniqueDel(self):
		mol = self.mc.molecule("enz3:mature", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.uniqueDel(newEnz3)
		self.assertEqual(mol.countUnique(), 0.0)

	@noseAttrib.attr('uniqueTest')
	def test_uniqueDel_objectNotCorrect(self):
		mol1 = self.mc.molecule("enz3:mature", "c")
		mol2 = self.mc.molecule("enz4:mature", "c")
		newEnz3 = mol1.uniqueNew()
		newEnz4 = mol2.uniqueNew()
		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
			mol1.uniqueDel(newEnz4)
		self.assertEqual(context.exception.message, 'Unique object to delete does not match row in unique table!\n')

	@noseAttrib.attr('uniqueTest')
	def test_uniquesWithAttrs(self):
		mol = self.mc.molecule("enz3:mature", "c")
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
		mol = self.mc.molecule("enz3:mature", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "A", "attr3" : "C"})
		with self.assertRaises(wcMoleculeCounts.uniqueException) as context:
			mol.uniquesWithAttrs({"attr1" : "A", "attr4" : "B"})
		self.assertEqual(context.exception.message, 'A specified attribute is not included in knoweldge base for this unique object!\n')

	@noseAttrib.attr('uniqueTest')
	def test_makeGetter(self):
		mol = self.mc.molecule("enz3:mature", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		self.assertEqual(newEnz3.attr1(), "A")
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")
		self.assertEqual(id(newEnz3), id(newEnz3.objects()))

	@noseAttrib.attr('uniqueTest')
	def test_makeSetter(self):
		mol = self.mc.molecule("enz3:mature", "c")
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		newEnz3.attr1Is("Z")
		self.assertEqual(newEnz3.attr1(), "Z")
		self.assertEqual(newEnz3.attr2(), "B")
		self.assertEqual(newEnz3.attr3(), "C")
		self.assertEqual(id(newEnz3), id(newEnz3.objects()))

	@noseAttrib.attr('uniqueTest')
	def test_dMassIs(self):
		mol = self.mc.molecule("enz3:mature", "c")
		mol.countBulkInc(1)
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.dMassIs(2)
		self.assertEqual(mol.massAll(), 8.0)

	@noseAttrib.attr('uniqueTest')
	def test_dMassInc(self):
		mol = self.mc.molecule("enz3:mature", "c")
		mol.countBulkInc(1)
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.dMassInc(2)
		self.assertEqual(mol.massAll(), 8.0)
		mol.dMassInc(1)
		self.assertEqual(mol.massAll(), 9.0)

	@noseAttrib.attr('uniqueTest')
	def test_dMassDec(self):
		mol = self.mc.molecule("enz3:mature", "c")
		mol.countBulkInc(1)
		newEnz3 = mol.uniqueNew({"attr1" : "A", "attr2" : "B", "attr3" : "C"})
		mol.dMassIs(2)
		self.assertEqual(mol.massAll(), 8.0)
		mol.dMassDec(1)
		self.assertEqual(mol.massAll(), 7.0)

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_sameAsDefault(self):
		mol = self.mc.molecule("enz4:mature", "c")
		e4 = mol.uniqueNew({"attr4_1" : "A", "attr4_2" : "B"})

		self.assertTrue(hasattr(e4, 'attr4_1'))
		self.assertTrue(hasattr(e4, 'attr4_1Is'))
		self.assertTrue(hasattr(e4, 'attr4_2'))
		self.assertTrue(hasattr(e4, 'attr4_2Is'))
		self.assertTrue(hasattr(e4, 'objects'))
		self.assertTrue(hasattr(e4, 'objectsIs'))

		self.assertEqual(e4.attr4_1(), "A")
		self.assertEqual(e4.attr4_2(), "B")
		self.assertEqual(id(e4.objects()), id(e4))
		e4.attr4_1Is("X")
		self.assertEqual(e4.attr4_1(), "X")
		e4.attr4_2Is("Y")
		self.assertEqual(e4.attr4_2(), "Y")

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_addNewFunction(self):
		mol = self.mc.molecule("enz5:mature", "c")
		e5 = mol.uniqueNew({"attr5_1" : "A"})
		self.assertEqual(e5.test_function(), "Created new function")

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_keepsDefaultFunctions(self):
		mol = self.mc.molecule("enz5:mature", "c")

		e5 = mol.uniqueNew({"attr5_1" : "A"})

		self.assertTrue(hasattr(e5, 'attr5_1'))
		self.assertTrue(hasattr(e5, 'attr5_1Is'))
		self.assertTrue(hasattr(e5, 'objects'))
		self.assertTrue(hasattr(e5, 'objectsIs'))

	@noseAttrib.attr('uniqueTest')
	def test_metaClass_multipleUniqueMetaclasses(self):
		mol4 = self.mc.molecule("enz4:mature", "c")
		e4 = mol4.uniqueNew({"attr4_1" : "A", "attr4_2" : "B"})
		mol5 = self.mc.molecule("enz5:mature", "c")
		e5 = mol5.uniqueNew({"attr5_1" : "A"})

		self.assertTrue(hasattr(e4, 'attr4_1'))
		self.assertTrue(hasattr(e4, 'attr4_1Is'))
		self.assertTrue(hasattr(e4, 'attr4_2'))
		self.assertTrue(hasattr(e4, 'attr4_2Is'))
		self.assertTrue(hasattr(e4, 'objects'))
		self.assertTrue(hasattr(e4, 'objectsIs'))

		self.assertTrue(hasattr(e5, 'attr5_1'))
		self.assertTrue(hasattr(e5, 'attr5_1Is'))
		self.assertTrue(hasattr(e5, 'objects'))
		self.assertTrue(hasattr(e5, 'objectsIs'))
		self.assertTrue(hasattr(e5, 'test_function'))


class enz4Unique_same(object):
	registrationId = "enz4:mature"
	__metaclass__ = wcMoleculeCounts.MoleculeUniqueMeta

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
	registrationId = "enz5:mature"
	__metaclass__ = wcMoleculeCounts.MoleculeUniqueMeta

	def test_function(self):
		return "Created new function"

