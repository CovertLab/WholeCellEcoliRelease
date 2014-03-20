'''
test_unique_objects_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.containers.unique_objects_container import UniqueObjectsContainer, UniqueObjectsContainerException

TEST_KB = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		},
	'DNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		}
	}


class Test_UniqueObjectsContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = createContainer()


	def tearDown(self):
		pass

	# Interface tests

	# Adding/removing objects
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_add_molecule(self):
		self.container.objectNew('RNA polymerase')

		self.assertEqual(len(self.container.objectsWithName_newMethod('RNA polymerase')), 21)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_add_molecules(self):
		self.container.objectsNew('RNA polymerase', 20)

		self.assertEqual(
			len(self.container.objectsWithName_newMethod('RNA polymerase')),
			40
			)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_delete_molecules(self):
		molecules = self.container.objectsWithName_newMethod('RNA polymerase')

		self.container.objectsDel(molecules)

		# TODO: reimplement
		# self.assertEqual(
		# 	self.container.objectsWithName_newMethod('RNA polymerase'),
		# 	set()
		# 	)

		# Make sure access to deleted molecules is blocked
		for molecule in molecules:
			# Weird logic to get one molecule
			break

		with self.assertRaises(UniqueObjectsContainerException) as context:
			molecule.attr('boundToChromosome')

		self.assertEqual(str(context.exception), 'Attempted to access an inactive molecule.')

		with self.assertRaises(UniqueObjectsContainerException) as context:
			molecule.attrIs(boundToChromosome = False)

		self.assertEqual(str(context.exception), 'Attempted to access an inactive molecule.')

	# Querying
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_empty_query(self):
		molecules = self.container.objectsWithName_newMethod('RNA polymerase')

		self.assertEqual(len(molecules), 20)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_bool_query(self):
		molecules = self.container.objectsWithName_newMethod(
			'RNA polymerase',
			boundToChromosome = ('==', False)
			)

		self.assertEqual(len(molecules), 10)

		for molecule in molecules:
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				False
				)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_numeric_query(self):
		molecules = self.container.objectsWithName_newMethod(
			'RNA polymerase',
			chromosomeLocation = ('>', 0)
			)

		self.assertEqual(len(molecules), 5)

		for molecule in molecules:
			self.assertGreater(
				molecule.attr('chromosomeLocation'),
				0
				)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_compound_query(self):
		molecules = self.container.objectsWithName_newMethod(
			'RNA polymerase',
			boundToChromosome = ('!=', False),
			chromosomeLocation = ('>', 0)
			)

		self.assertEqual(len(molecules), 5)

		for molecule in molecules:
			self.assertGreater(
				molecule.attr('chromosomeLocation'),
				0
				)

	# Attribute access
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_attribute_setting(self):
		for molecule in self.container.objectsWithName_newMethod('RNA polymerase'):
			molecule.attrIs(
				boundToChromosome = True,
				chromosomeLocation = 100
				)

		for molecule in self.container.objectsWithName_newMethod('RNA polymerase'):
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

			self.assertEqual(
				molecule.attr('chromosomeLocation'),
				100
				)

	# Internal tests

	# Bookkeeping attributes
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_time_setting(self):
		self.container._timeIs(50)

		allMolecules = self.container.objectsWithName_newMethod('RNA polymerase')
		newTime = self.container.objectsWithName_newMethod('RNA polymerase', _time = ('==', 50))

		# TODO: reimplement
		#self.assertEqual(allMolecules, newTime)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_deleted_entry_flushing(self):
		# First, make sure that deleted entries are not overwritten
		molecules = self.container.objectsWithName_newMethod('RNA polymerase')
		indexes = {molecule._objectIndex for molecule in molecules}

		self.container.objectsDel(molecules)

		newMolecule = self.container.objectNew(
			'RNA polymerase'
			)

		self.assertTrue(newMolecule._objectIndex not in indexes)

		# Next, flush the deleted entries and confirm that an old entry is overwritten
		self.container._flushDeleted()
		newMolecule = self.container.objectNew(
			'RNA polymerase'
			)

		self.assertTrue(newMolecule._objectIndex in indexes)

	# Global references
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_global_index_mapping(self):
		globalArray = self.container._arrays[self.container._globalRefIndex]

		for molecule in self.container.objectsWithName_newMethod('RNA polymerase'):
			globalIndex = molecule.attr('_globalIndex')

			globalEntry = globalArray[globalIndex]

			arrayIndex = globalEntry['_arrayIndex']
			objectIndex = globalEntry['_objectIndex']

			self.assertEqual(molecule._arrayIndex, arrayIndex)
			self.assertEqual(molecule._objectIndex, objectIndex)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_global_index_removal(self):
		globalArray = self.container._arrays[self.container._globalRefIndex]

		for molecule in self.container.objectsWithName_newMethod('RNA polymerase'):
			# Weird logic to get one molecule
			break

		globalIndex = molecule.attr('_globalIndex')
		
		self.container.objectDel(molecule)
		self.container._flushDeleted()

		globalEntry = globalArray[globalIndex]

		deletedEntry = np.zeros(1, dtype = globalArray.dtype)

		self.assertEqual(globalEntry, deletedEntry)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_eq_method(self):
		# Test against self
		self.assertEqual(self.container, self.container)

		# Test against other, presumably identical container
		otherContainer = createContainer()

		self.assertEqual(self.container, otherContainer)


	# Testing new methods TODO: migrate

	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objects(self):
		objectSet = self.container.objects_newMethod()

		self.assertEqual(len(objectSet), 20)

		for obj in objectSet:
			self.assertIn(obj, objectSet)

		objectSet = self.container.objects_newMethod(chromosomeLocation = ('>', 0))

		self.assertEqual(len(objectSet), 5)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectsWithName(self):
		self.container.objectsNew('DNA polymerase', 20)

		objectSet = self.container.objectsWithName_newMethod('RNA polymerase')

		self.assertEqual(len(objectSet), 20)

		for obj in objectSet:
			self.assertIn(obj, objectSet)

		objectSet = self.container.objectsWithName_newMethod(
			'RNA polymerase', chromosomeLocation = ('>', 0))

		self.assertEqual(len(objectSet), 5)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectsWithNames(self):
		self.container.objectsNew('DNA polymerase', 20)

		objectSet = self.container.objectsWithNames_newMethod(
			['RNA polymerase', 'DNA polymerase'])

		self.assertEqual(len(objectSet), 40)

		for obj in objectSet:
			self.assertIn(obj, objectSet)

		objectSet = self.container.objectsWithNames_newMethod(
			['RNA polymerase', 'DNA polymerase'], chromosomeLocation = ('==', 0))

		self.assertEqual(len(objectSet), 35)


def createContainer():
	container = UniqueObjectsContainer(TEST_KB)
	
	container.objectsNew(
		'RNA polymerase',
		10,
		)

	container.objectsNew(
		'RNA polymerase',
		5,
		boundToChromosome = True,
		chromosomeLocation = 0
		)

	container.objectsNew(
		'RNA polymerase',
		5,
		boundToChromosome = True,
		chromosomeLocation = 50
		)

	return container
