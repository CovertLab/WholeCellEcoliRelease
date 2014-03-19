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

import wholecell.utils.unique_objects_container

TEST_KB = {
	'RNA polymerase':{
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

		self.assertEqual(len(self.container.objects('RNA polymerase')), 21)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_add_molecules(self):
		self.container.objectsNew('RNA polymerase', 20)

		self.assertEqual(
			len(self.container.objects('RNA polymerase')),
			40
			)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_delete_molecules(self):
		molecules = self.container.objects('RNA polymerase')

		self.container.objectsDel(molecules)

		self.assertEqual(
			self.container.objects('RNA polymerase'),
			set()
			)

		# Make sure access to deleted molecules is blocked
		molecule = molecules.pop()

		with self.assertRaises(Exception) as context:
			molecule.attr('boundToChromosome')

		self.assertEqual(str(context.exception), 'Attempted to access an inactive molecule.')

		with self.assertRaises(Exception) as context:
			molecule.attrIs('boundToChromosome', False)

		self.assertEqual(str(context.exception), 'Attempted to access an inactive molecule.')

	# Querying
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_empty_query(self):
		molecules = self.container.evaluateQuery('RNA polymerase')

		self.assertEqual(len(molecules), 20)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_bool_query(self):
		molecules = self.container.evaluateQuery(
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
		molecules = self.container.evaluateQuery(
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
		molecules = self.container.evaluateQuery(
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
		for molecule in self.container.iterObjects('RNA polymerase'):
			molecule.attrIs('boundToChromosome', True)
			molecule.attrIs('chromosomeLocation', 100)

		for molecule in self.container.iterObjects('RNA polymerase'):
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

			self.assertEqual(
				molecule.attr('chromosomeLocation'),
				100
				)

	# Query objects
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_query_objects(self):
		query = self.container.queryNew('RNA polymerase', boundToChromosome = ('==', True))

		self.container.updateQueries()

		for molecule in query.iterObjects():
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

		for molecule in query.iterObjects():
			molecule.attrIs('boundToChromosome', False)

		self.container.updateQueries()

		self.assertEqual(query.objects(), set())

	# Set operations
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_molecule_set_operations(self):
		allMolecules = self.container.objects('RNA polymerase')

		chromosomeBound = self.container.evaluateQuery(
			'RNA polymerase',
			boundToChromosome = ('==', True)
			)

		self.assertEqual(allMolecules & allMolecules, allMolecules)

		self.assertEqual(allMolecules & chromosomeBound, chromosomeBound)

		self.assertEqual(allMolecules | chromosomeBound, allMolecules)

		self.assertEqual(len(allMolecules - chromosomeBound), 10)

	# Internal tests

	# Bookkeeping attributes
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_time_setting(self):
		self.container._timeIs(50)

		allMolecules = self.container.objects('RNA polymerase')
		newTime = self.container.evaluateQuery('RNA polymerase', _time = ('==', 50))

		self.assertEqual(allMolecules, newTime)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_deleted_entry_flushing(self):
		# First, make sure that deleted entries are not overwritten
		molecules = self.container.objects('RNA polymerase')
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

		for molecule in self.container.objects('RNA polymerase'):
			globalIndex = molecule.attr('_globalIndex')

			globalEntry = globalArray[globalIndex]

			arrayIndex = globalEntry['_arrayIndex']
			objectIndex = globalEntry['_objectIndex']

			self.assertEqual(molecule._arrayIndex, arrayIndex)
			self.assertEqual(molecule._objectIndex, objectIndex)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_global_index_removal(self):
		globalArray = self.container._arrays[self.container._globalRefIndex]

		molecule = self.container.objects('RNA polymerase').pop()

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


def createContainer():
	container = wholecell.utils.unique_objects_container.UniqueObjectsContainer(
		TEST_KB)
	
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
