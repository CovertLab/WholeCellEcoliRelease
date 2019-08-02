'''
test_unique_objects_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import absolute_import, division, print_function

import cPickle
import unittest
import os
import shutil
import tempfile

import numpy as np
import numpy.testing as npt
import nose.plugins.attrib as noseAttrib

from wholecell.containers.unique_objects_container import (
	UniqueObjectsContainer, UniqueObjectsContainerException,
	UniqueObjectsPermissionException, Access
	)

from wholecell.io.tablereader import TableReader
from wholecell.io.tablewriter import TableWriter


TEST_KB = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		},
	'DNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		},
	'Chocolate':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32',
		'percent':'float32',
		'nuts':'bool',
		},
	}


class Test_UniqueObjectsContainer(unittest.TestCase):
	def setUp(self):
		self.container = createContainer()
		self.test_dir = None

	def tearDown(self):
		if self.test_dir:
			shutil.rmtree(self.test_dir)

	def make_test_dir(self):
		if not self.test_dir:
			self.test_dir = tempfile.mkdtemp()


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
		molecules = self.container.objectsInCollection(
			'RNA polymerase', access=Access.READ_EDIT_DELETE)

		self.container.objectsDel(molecules)

		self.assertEqual(
			len(self.container.objectsInCollection('RNA polymerase')),
			0
			)

		# Make sure access to deleted molecules is blocked
		for molecule in molecules:
			# Weird logic to get one molecule
			break
		else:
			raise ValueError("Expected at least one molecule")

		with self.assertRaises(UniqueObjectsContainerException) as context:
			molecule.attr('boundToChromosome')

		self.assertEqual(str(context.exception), 'Attempted to access an inactive object.')

		with self.assertRaises(UniqueObjectsContainerException) as context:
			molecule.attrIs(boundToChromosome = False)

		self.assertEqual(str(context.exception), 'Attempted to access an inactive object.')


	# TODO(jerry): Test _UniqueObjectSet.delByIndexes().


	# Querying
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_query(self):
		molecules = self.container.objectsInCollection('RNA polymerase')

		self.assertEqual(len(molecules), 20)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_merge(self):
		container = createContainer()
		n_objects = len(container.objects())

		container.add_request(
			type="new_molecule",
			collectionName="RNA polymerase",
			nObjects=1,
			attributes={
				"boundToChromosome": True,
				"chromosomeLocation": 0
				}
			)

		# Object should not be added before .merge()
		self.assertEqual(n_objects, len(container.objects()))

		container.merge()
		self.assertEqual(n_objects + 1, len(container.objects()))


	# Attribute access
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_attribute_setting(self):
		for molecule in self.container.objectsInCollection('RNA polymerase', access=Access.READ_EDIT):
			molecule.attrIs(
				boundToChromosome = True,
				chromosomeLocation = 100,
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


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectsInCollection(self):
		self.container.objectsNew('DNA polymerase', 20)

		objectSet = self.container.objectsInCollection('RNA polymerase')

		self.assertEqual(len(objectSet), 20)

		for obj in objectSet:
			self.assertIn(obj, objectSet)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectsInCollections(self):
		self.container.objectsNew('DNA polymerase', 20)

		objectSet = self.container.objectsInCollections(
			['RNA polymerase', 'DNA polymerase'])

		self.assertEqual(len(objectSet), 40)

		for obj in objectSet:
			self.assertIn(obj, objectSet)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_counts(self):
		self.container.objectsNew('DNA polymerase', 15)

		npt.assert_array_equal((0, 15, 20), self.container.counts(), str(self.container.objectNames()))
		self.assertEqual(15, self.container.counts(["DNA polymerase"]), "DNA polymerase")

		self.container.objectsNew('Chocolate', 1000, nuts=True, percent=95.0)
		npt.assert_array_equal((1000, 15, 20), self.container.counts(), str(self.container.objectNames()))
		npt.assert_array_equal((20, 15), self.container.counts(["RNA polymerase", "DNA polymerase"]), "RNA, DNA")


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
		else:
			raise ValueError("expected at least one molecule")

		globalIndex = molecule.attr('_globalIndex')

		self.container.objectDel(molecule)

		globalEntry = globalArray[globalIndex]

		deletedEntry = np.zeros(1, dtype = globalArray.dtype)

		self.assertEqual(globalEntry, deletedEntry)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_unique_index(self):
		objectSet = self.container.objects()
		unique_index = objectSet.attr("_uniqueIndex")

		# Check that all unique index values are unique
		self.assertEqual(len(np.unique(unique_index)), len(unique_index))

		# Add new objects to container
		self.container.objectsNew('Chocolate', 10, nuts=True, percent=95.0)
		objectSet = self.container.objects()
		unique_index = objectSet.attr("_uniqueIndex")

		# Check that all unique index values are unique
		self.assertEqual(len(np.unique(unique_index)), len(unique_index))


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_eq_method(self):
		# Test against self
		self.assertEqual(self.container, self.container)
		self.assertNotEqual(self.container, object())

		# Test against other container with different massdiff names
		container_empty_massdiff = createContainer()
		self.assertEqual(self.container, container_empty_massdiff)

		container_empty_massdiff.submass_diff_names_list = []
		self.assertNotEqual(self.container, container_empty_massdiff)

		# Test against other, presumably identical container
		otherContainer = createContainer()
		self.assertEqual(self.container, otherContainer)

		self.container.objectsNew(
			'Chocolate', 3, percent=90.0, nuts=True, chromosomeLocation=2001)
		self.assertNotEqual(self.container, otherContainer)

		otherContainer.objectsNew(
			'Chocolate', 3, chromosomeLocation=2001, percent=90.0, nuts=True)
		self.assertEqual(self.container, otherContainer)
		npt.assert_array_equal(self.container._globalReference, otherContainer._globalReference)

		self.container.objectsNew(
			'Chocolate', 1, percent=80.0, chromosomeLocation=2001, nuts=False)
		self.assertNotEqual(self.container, otherContainer)

		otherContainer.objectsNew(
			'Chocolate', 1, nuts=False, percent=80.0, chromosomeLocation=2001)
		self.assertEqual(self.container, otherContainer)

		self.container.objectsNew(
			'Chocolate', 1, percent=95.0, chromosomeLocation=3000, nuts=False)
		self.assertNotEqual(self.container, otherContainer)

		otherContainer.objectsNew(
			'Chocolate', 1, nuts=False, percent=95.0, chromosomeLocation=1000)
		self.assertNotEqual(self.container, otherContainer)


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_emptyLike(self):
		e1 = self.container.emptyLike()
		e2 = self.container.emptyLike()
		self.assertNotEqual(self.container, e1)
		self.assertEqual(e1, e2)
		self.assertEqual(tuple(sorted(TEST_KB.keys())), e1.objectNames())
		self.assertEqual(e1._specifications, e2._specifications)
		npt.assert_array_equal((0, 0, 0), e1.counts())

		self.assertEqual(0,
			np.where(e1._globalReference['_entryState'] == e1._entryActive)[0].size)
		self.assertEqual(0, len([obj for obj in e1.objects()]))
		self.assertEqual(0, len([obj for obj in e1.objectsInCollection('Chocolate')]))
		self.assertEqual(0, len([obj for obj in e1.objectsInCollections(TEST_KB.keys())]))

	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_copyContents(self):
		self.container.objectsNew('Chocolate', 17, percent=88.8)
		npt.assert_array_equal((17, 0, 20), self.container.counts())

		e1 = self.container.emptyLike()
		self.assertNotEqual(self.container, e1)
		e1.loadSnapshot(self.container)
		self.assertEqual(self.container, e1)
		npt.assert_array_equal((17, 0, 20), e1.counts())

		KB2 = {'RNA': {'boundToChromosome': 'bool', 'chromosomeLocation': 'uint32'},
			'RSA:': {'boundToChromosome': 'bool', 'chromosomeLocation': 'uint32'}}
		c1 = UniqueObjectsContainer(KB2, [])
		self.assertNotEqual(self.container, c1)

		with self.assertRaises(ValueError):  # different specifications
			c1.loadSnapshot(self.container)


	# I/O

	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_pickle(self):
		data = cPickle.dumps(self.container)
		container2 = cPickle.loads(data)
		self.assertEqual(self.container, container2)

		# print("Pickled a UniqueObjectsContainer to {} bytes".format(len(data)))

		# Test that internal fields are properly restored which __eq__() assumes.
		self.assertEqual(self.container.objectNames(), container2.objectNames())
		self.assertEqual(self.container._nameToIndexMapping, container2._nameToIndexMapping)
		self.assertEqual(self.container.submass_diff_names_list, container2.submass_diff_names_list)
		npt.assert_array_equal(self.container._globalReference, container2._globalReference)

		# Test on a collection with added & deleted entries.
		self.container.objectsNew('Chocolate', 101, percent=95.0, nuts=False)
		self.assertNotEqual(self.container, container2)

		# Test that container with unapplied requests does not get pickled
		self.container.add_request(
				type="delete",
				globalIndexes=np.array([0]),
			)
		with self.assertRaises(UniqueObjectsContainerException) as context:
			cPickle.dumps(self.container)

		self.assertEqual(
			str(context.exception),
			"Cannot pickle container with unapplied requests. Run .merge() to apply the requests before pickling."
			)

		self.container._requests = []

		counts = self.container.counts(['RNA polymerase', 'DNA polymerase'])
		rna_objects = self.container.objectsInCollection('RNA polymerase')
		self.container.objectDel(rna_objects[0])
		npt.assert_array_equal(counts - [1, 0],
			self.container.counts(['RNA polymerase', 'DNA polymerase']))

		data = cPickle.dumps(self.container)
		container3 = cPickle.loads(data)
		self.assertEqual(self.container, container3)
		self.assertEqual(self.container.objectNames(), container3.objectNames())
		self.assertEqual(self.container._nameToIndexMapping, container3._nameToIndexMapping)
		npt.assert_array_equal(self.container._globalReference, container3._globalReference)

		# print("Pickled a UniqueObjectsContainer to {} bytes".format(len(data)))

		# Test that an unpickled container is mutable as always. It's at risk
		# because np.frombuffer(f.read()) gets an ndarray onto a readonly byte str.
		container3.objectsNew('Chocolate', 3, percent=75.0, nuts=True)
		self.assertEqual((104,), container3.counts(['Chocolate']))


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_write_table(self):
		"""Test writing the container to a Table."""
		self.make_test_dir()
		path = os.path.join(self.test_dir, 'UniqueObjects')

		self.container.objectsNew('Chocolate', 2001, percent=75.0, nuts=True)

		rna_objects = self.container.objectsInCollection('RNA polymerase')
		self.assertEqual(20, len(rna_objects))
		self.container.objectDel(rna_objects[1])
		rna_objects2 = self.container.objectsInCollection('RNA polymerase')
		self.assertEqual(19, len(rna_objects2))

		table_writer = TableWriter(path)
		self.container.tableCreate(table_writer)
		self.container.tableAppend(table_writer)
		table_writer.close()

		# Read and check the data.
		# ASSUMES: One table row appended, readColumn() squeezes the array to
		# 1D, and various internal details about how this container writes a
		# table such as including unused array elements.
		table_reader = TableReader(path)
		names = self.container.objectNames()
		self.assertEqual(list(names), table_reader.readAttribute('collectionNames'))
		for index, name in enumerate(names):
			expected = self.container._collections[index]
			npt.assert_array_equal(expected, table_reader.readColumn(name))


	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_objectSet_attribute_accessing(self):
		objectSet = self.container.objects(
			access=Access.READ_EDIT
			)

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

	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_read_only(self):
		molecules = self.container.objectsInCollection(
			'RNA polymerase', access=Access.READ_ONLY
			)
		n_molecules = len(molecules)

		with self.assertRaises(UniqueObjectsPermissionException) as context:
			molecules.attrIs(chromosomeLocation = np.ones(n_molecules))

		self.assertEqual(
			str(context.exception),
			"Can't modify attributes of read-only objects."
			)

		with self.assertRaises(UniqueObjectsPermissionException) as context:
			molecules.delByIndexes(np.ones(n_molecules))

		self.assertEqual(
			str(context.exception),
			"Can't delete molecules from read-only or read-and-edit only objects."
			)

	@noseAttrib.attr('smalltest', 'uniqueObjects', 'containerObject')
	def test_read_edit(self):
		molecules = self.container.objectsInCollection(
			'RNA polymerase', access=Access.READ_EDIT
			)
		n_molecules = len(molecules)

		with self.assertRaises(UniqueObjectsPermissionException) as context:
			molecules.delByIndexes(np.ones(n_molecules))

		self.assertEqual(
			str(context.exception),
			"Can't delete molecules from read-only or read-and-edit only objects."
			)



def createContainer():
	container = UniqueObjectsContainer(TEST_KB, ["massDiff_RNA"])

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
