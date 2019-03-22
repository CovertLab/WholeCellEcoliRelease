'''
test_unique_objects_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import division

import unittest
import os

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

		self.assertEqual(len(self.container.objectsInCollection('RNA polymerase')), 21)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_add_molecules(self):
		self.container.objectsNew('RNA polymerase', 20)

		self.assertEqual(
			len(self.container.objectsInCollection('RNA polymerase')),
			40
			)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_delete_molecules(self):
		molecules = self.container.objectsInCollection('RNA polymerase')

		self.container.objectsDel(molecules)

		self.assertEqual(
			len(self.container.objectsInCollection('RNA polymerase')),
			0
			)

		# Make sure access to deleted molecules is blocked
		for molecule in molecules:
			# Weird logic to get one molecule
			break

		with self.assertRaises(UniqueObjectsContainerException) as context:
			molecule.attr('boundToChromosome')

		self.assertEqual(str(context.exception), 'Attempted to access an inactive object.')

		with self.assertRaises(UniqueObjectsContainerException) as context:
			molecule.attrIs(boundToChromosome = False)

		self.assertEqual(str(context.exception), 'Attempted to access an inactive object.')

	# Querying
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_empty_query(self):
		molecules = self.container.objectsInCollection('RNA polymerase')

		self.assertEqual(len(molecules), 20)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_bool_query(self):
		molecules = self.container.objectsInCollection(
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
		molecules = self.container.objectsInCollection(
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
		molecules = self.container.objectsInCollection(
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


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_contains_query(self):
		validLocations = np.array([0, 50])

		molecules = self.container.objectsInCollection(
			'RNA polymerase',
			boundToChromosome = ("==", True),
			chromosomeLocation = ('in', validLocations)
			)

		self.assertEqual(len(molecules), 10)

		for molecule in molecules:
			self.assertIn(
				molecule.attr('chromosomeLocation'),
				validLocations
				)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_not_contains_query(self):
		invalidLocations = np.array([0])

		molecules = self.container.objectsInCollection(
			'RNA polymerase',
			boundToChromosome = ("==", True),
			chromosomeLocation = ('not in', invalidLocations)
			)

		self.assertEqual(len(molecules), 5)

		for molecule in molecules:
			self.assertNotIn(
				molecule.attr('chromosomeLocation'),
				invalidLocations
				)

	# Attribute access
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_attribute_setting(self):
		for molecule in self.container.objectsInCollection('RNA polymerase'):
			molecule.attrIs(
				boundToChromosome = True,
				chromosomeLocation = 100
				)

		for molecule in self.container.objectsInCollection('RNA polymerase'):
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

			self.assertEqual(
				molecule.attr('chromosomeLocation'),
				100
				)

	# Object access

	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objects(self):
		objectSet = self.container.objects()

		self.assertEqual(len(objectSet), 20)

		for obj in objectSet:
			self.assertIn(obj, objectSet)

		objectSet = self.container.objects(chromosomeLocation = ('>', 0))

		self.assertEqual(len(objectSet), 5)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectsInCollection(self):
		self.container.objectsNew('DNA polymerase', 20)

		objectSet = self.container.objectsInCollection('RNA polymerase')

		self.assertEqual(len(objectSet), 20)

		for obj in objectSet:
			self.assertIn(obj, objectSet)

		objectSet = self.container.objectsInCollection(
			'RNA polymerase', chromosomeLocation = ('>', 0))

		self.assertEqual(len(objectSet), 5)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectsInCollections(self):
		self.container.objectsNew('DNA polymerase', 20)

		objectSet = self.container.objectsInCollections(
			['RNA polymerase', 'DNA polymerase'])

		self.assertEqual(len(objectSet), 40)

		for obj in objectSet:
			self.assertIn(obj, objectSet)

		objectSet = self.container.objectsInCollections(
			['RNA polymerase', 'DNA polymerase'], chromosomeLocation = ('==', 0))

		self.assertEqual(len(objectSet), 35)

	# Internal tests

	# Global references
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_global_index_mapping(self):
		globalArray = self.container._globalReference

		for molecule in self.container.objectsInCollection('RNA polymerase'):
			globalIndex = molecule.attr('_globalIndex')

			globalEntry = globalArray[globalIndex]

			arrayIndex = globalEntry['_collectionIndex']
			objectIndex = globalEntry['_objectIndex']

			self.assertEqual(molecule._collectionIndex, arrayIndex)
			self.assertEqual(molecule._objectIndex, objectIndex)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_global_index_removal(self):
		globalArray = self.container._globalReference

		for molecule in self.container.objectsInCollection('RNA polymerase'):
			# Weird logic to get one molecule
			break

		globalIndex = molecule.attr('_globalIndex')

		self.container.objectDel(molecule)

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

	# TODO: test saving/loading

	# @noseAttrib.attr('mediumtest', 'uniqueObjects', 'containerObject', 'saveload')
	# def test_save_load(self):
	# 	# Create file, save values, close
	# 	path = os.path.join('fixtures', 'test', 'test_unique_objects_container.hdf')

	# 	h5file = tables.open_file(
	# 		path,
	# 		mode = 'w',
	# 		title = 'File for UniqueObjectsContainer IO'
	# 		)

	# 	self.container.pytablesCreate(h5file)

	# 	self.container.timeStepIs(0)

	# 	self.container.pytablesAppend(h5file)

	# 	h5file.close()

	# 	# Open, load, and compare
	# 	h5file = tables.open_file(path)

	# 	loadedContainer = UniqueObjectsContainer(TEST_KB)

	# 	loadedContainer.pytablesLoad(h5file, 0)

	# 	self.assertEqual(
	# 		self.container,
	# 		loadedContainer
	# 		)

	# 	h5file.close()


	# @noseAttrib.attr('mediumtest', 'uniqueObjects', 'containerObject', 'saveload')
	# def test_save_load_later(self):
	# 	# Create file, save values, close
	# 	path = os.path.join('fixtures', 'test', 'test_unique_objects_container.hdf')

	# 	h5file = tables.open_file(
	# 		path,
	# 		mode = 'w',
	# 		title = 'File for UniqueObjectsContainer IO'
	# 		)

	# 	self.container.pytablesCreate(h5file)

	# 	self.container.timeStepIs(0)

	# 	self.container.pytablesAppend(h5file)

	# 	self.container.objectsNew('DNA polymerase', 5)

	# 	self.container.timeStepIs(1)

	# 	self.container.pytablesAppend(h5file)

	# 	h5file.close()

	# 	# Open, load, and compare
	# 	h5file = tables.open_file(path)

	# 	loadedContainer = UniqueObjectsContainer(TEST_KB)

	# 	loadedContainer.pytablesLoad(h5file, 1)

	# 	self.assertEqual(
	# 		self.container,
	# 		loadedContainer
	# 		)

	# 	h5file.close()


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectSet_attribute_accessing(self):
		objectSet = self.container.objects(chromosomeLocation = ('>=', 0))

		self.assertEqual(len(objectSet), 20)

		# Test getters

		chromosomeLocation = objectSet.attr('chromosomeLocation')

		self.assertEqual(
			{0, 50},
			set(chromosomeLocation)
			)

		self.assertEqual(
			15,
			(chromosomeLocation == 0).sum()
			)

		self.assertEqual(
			5,
			(chromosomeLocation == 50).sum()
			)

		chromosomeLocation, boundToChromosome = objectSet.attrs(
			'chromosomeLocation', 'boundToChromosome')

		self.assertEqual(
			10,
			boundToChromosome.sum()
			)

		# Test setters

		objectSet.attrIs(chromosomeLocation = np.ones(20))

		chromosomeLocation = objectSet.attr('chromosomeLocation')

		self.assertEqual(
			20,
			(chromosomeLocation == 1).sum()
			)

		objectSet.attrIs(
			chromosomeLocation = np.zeros(20),
			boundToChromosome = np.zeros(20, np.bool)
			)

		chromosomeLocation, boundToChromosome = objectSet.attrs(
			'chromosomeLocation', 'boundToChromosome')

		self.assertEqual(
			20,
			(chromosomeLocation == 0).sum()
			)

		self.assertEqual(
			20,
			(~boundToChromosome).sum()
			)


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
