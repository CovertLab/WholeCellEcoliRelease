'''
test_bulk_objects_container.py

@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os
import shutil
import tempfile
from typing import Iterable
import unittest

import numpy as np
import numpy.testing as npt

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.io.tablereader import TableReader
from wholecell.io.tablewriter import TableWriter


OBJECT_NAMES = ('ATP', 'glucose', 'glycine')
OBJECT_COUNTS = [100, 20, 10]

# A compressible object name list.
ELEMENTS = '''Actinium Aluminum Americium Barium Berkelium Beryllium
	Bohrium Cadmium Calcium Californium Cerium Cesium Chromium
	Copernicium Curium Darmstadtium Dubnium Dysprosium Einsteinium
	Erbium Europium Fermium Flerovium Francium Gadolinium Gallium Germanium
	Hafnium Hassium Helium Holmium Indium Iridium Lanthanum Lawrencium
	Lithium Livermorium Lutetium Magnesium Meitnerium Mendelevium Molybdenum
	Moscovium Neodymium Neptunium Nihonium Niobium Nobelium Osmium Palladium
	Platinum Plutonium Polonium Potassium Praseodymium Promethium Protactinium
	Radium Rhenium Rhodium Roentgenium Rubidium Ruthenium Rutherfordium
	Samarium Scandium Seaborgium Selenium Sodium Strontium Tantalum Technetium
	Tellurium Terbium Thallium Thorium Thulium Titanium Uranium Vanadium
	Ytterbium Yttrium Zirconium'''.split()


class Test_BulkObjectsContainer(unittest.TestCase):
	def setUp(self):
		self.container = createContainer()
		self.test_dir = None

	def tearDown(self):
		if self.test_dir:
			shutil.rmtree(self.test_dir)

	def make_test_dir(self):
		if not self.test_dir:
			self.test_dir = tempfile.mkdtemp()


	# Interface methods

	def test_counts(self):
		npt.assert_equal(self.container.counts(), OBJECT_COUNTS)

		npt.assert_equal(
			self.container.counts(OBJECT_NAMES[1:]),
			OBJECT_COUNTS[1:]
			)


	def test_countsIs(self):
		newCounts = [10, 20.0, 30]
		self.container.countsIs(newCounts)

		npt.assert_equal(self.container.counts(), newCounts)

		newCounts = [15, 5.0]
		self.container.countsIs(newCounts, OBJECT_NAMES[1:])

		npt.assert_equal(self.container.counts(OBJECT_NAMES[1:]), newCounts)


	def test_countsInc(self):
		incCounts = [10, 20.0, 30]
		newCounts = self.container.counts() + incCounts
		self.container.countsInc(incCounts)

		npt.assert_equal(self.container.counts(), newCounts)

		incCounts = [15, 5.0]
		newCounts = self.container.counts(OBJECT_NAMES[1:]) + incCounts
		self.container.countsInc(incCounts, OBJECT_NAMES[1:])

		npt.assert_equal(self.container.counts(OBJECT_NAMES[1:]), newCounts)


	def test_countsDec(self):
		decCounts = [30, 10.0, 5]
		newCounts = np.array(OBJECT_COUNTS) - np.array(decCounts)
		self.container.countsDec(decCounts)

		npt.assert_equal(self.container.counts(), newCounts)

		decCounts = [5.0, 2]
		newCounts = self.container.counts(OBJECT_NAMES[1:]) - np.array(decCounts)
		self.container.countsDec(decCounts, OBJECT_NAMES[1:])

		npt.assert_equal(self.container.counts(OBJECT_NAMES[1:]), newCounts)


	def test_dtype_float32(self):
		"""A BulkObjectsContainer with dtype=np.float32 should support
		fractional counts and deltas.
		"""
		container = BulkObjectsContainer(OBJECT_NAMES, dtype=np.float32)
		initialCounts = [10, 10.5, 20]
		container.countsIs(initialCounts)
		npt.assert_equal(container.counts(), initialCounts)

		incCounts = [10, 20.5, 30.5]
		newCounts = [20, 31, 50.5]
		container.countsInc(incCounts)
		npt.assert_equal(container.counts(), newCounts)

		decCounts = [1.5, 2, 3.5]
		newCounts = [18.5, 29, 47]
		container.countsDec(decCounts)
		npt.assert_equal(container.counts(), newCounts)

		countsView = container.countsView()
		newCounts = [28.5, 49.5, 77.5]
		countsView.countsInc(incCounts)
		npt.assert_equal(countsView.counts(), newCounts)

		newCounts = [27, 47.5, 74]
		countsView.countsDec(decCounts)
		npt.assert_equal(countsView.counts(), newCounts)


	def test_countsView_countsInc(self):
		countsView = self.container.countsView()
		incCounts = [10, 20.0, 30]
		newCounts = self.container.counts() + np.array(incCounts)
		countsView.countsInc(incCounts)

		npt.assert_equal(countsView.counts(), newCounts)

		countsView = self.container.countsView(OBJECT_NAMES[1:])
		incCounts = [15, 5.0]
		newCounts = self.container.counts(OBJECT_NAMES[1:]) + np.array(incCounts)
		countsView.countsInc(incCounts)

		npt.assert_equal(countsView.counts(), newCounts)


	def test_countsView_countsDec(self):
		countsView = self.container.countsView()
		decCounts = [30, 10.0, 5]
		newCounts = np.array(OBJECT_COUNTS) - np.array(decCounts)
		countsView.countsDec(decCounts)

		npt.assert_equal(countsView.counts(), newCounts)

		countsView = self.container.countsView(OBJECT_NAMES[1:])
		decCounts = [5.0, 2]
		newCounts = self.container.counts(OBJECT_NAMES[1:]) - np.array(decCounts)
		countsView.countsDec(decCounts)

		npt.assert_equal(countsView.counts(), newCounts)


	def test_count(self):
		npt.assert_equal(
			self.container.count(OBJECT_NAMES[0]),
			OBJECT_COUNTS[0]
			)


	def test_countIs(self):
		newCount = 40
		self.container.countIs(newCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)


	def test_countInc(self):
		incCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)

		incCount = 40.0
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)


	def test_countDec(self):
		decCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)

		decCount = 40.0
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)


	def test_emptyLike(self):
		e1 = self.container.emptyLike()
		e2 = self.container.emptyLike()
		self.assertNotEqual(self.container, e1)
		self.assertEqual(e1, e2)
		self.assertEqual(OBJECT_NAMES, e1.objectNames())
		npt.assert_array_equal((0, 0, 0), e1.counts())

	def test_copyCounts(self):
		npt.assert_array_equal(OBJECT_COUNTS, self.container.counts())

		c1 = self.container.emptyLike()
		c2 = BulkObjectsContainer(OBJECT_NAMES[1:])
		self.assertNotEqual(self.container, c1)
		self.assertNotEqual(c1, c2)

		c1.loadSnapshot(self.container)
		self.assertEqual(self.container, c1)
		npt.assert_array_equal(OBJECT_COUNTS, c1.counts())

		with self.assertRaises(ValueError):  # different object names
			c2.loadSnapshot(self.container)


	# Internal methods

	def test_namesToIndexes(self):
		# Test normal ordering
		names = OBJECT_NAMES
		npt.assert_equal(
			self.container._namesToIndexes(names),
			np.arange(len(OBJECT_NAMES))
			)

		# Test reverse ordering
		names = OBJECT_NAMES[::-1]
		npt.assert_equal(
			self.container._namesToIndexes(names),
			np.arange(len(OBJECT_NAMES))[::-1]
			)

		# Test subset
		names = OBJECT_NAMES[::2]
		npt.assert_equal(
			self.container._namesToIndexes(names),
			np.arange(len(OBJECT_NAMES))[::2]
			)


	def test_eq(self):
		self.assertEqual(self.container, self.container)
		self.assertNotEqual(self.container, object())

		newContainer = createContainer()

		self.assertEqual(newContainer, self.container)
		npt.assert_equal(self.container.counts(), newContainer.counts())

		container2 = BulkObjectsContainer(tuple('abcdefghijklmnopqrstuvwxyz'))
		assert self.container != container2  # different names and shape

		# __eq__() tests the counts. The dtypes may differ.
		newContainer.countsIs(7 + 9876 * np.arange(len(OBJECT_NAMES)))
		container3 = BulkObjectsContainer(OBJECT_NAMES, 'int16')
		self.assertNotEqual(newContainer, container3)
		container3.countsIs(7 + 9876 * np.arange(len(OBJECT_NAMES)))
		self.assertEqual(newContainer, container3)


	# I/O

	def test_pickle(self):
		# Test a float64 container.
		container = BulkObjectsContainer(ELEMENTS, dtype=np.float64)
		container.countsIs(1 + 1234.56 * np.arange(len(ELEMENTS)))

		data = cPickle.dumps(container, cPickle.HIGHEST_PROTOCOL)
		container2 = cPickle.loads(data)
		# print("Pickled a BulkObjectsContainer of {} float64 to {} bytes".format(
		# 	len(ELEMENT_NAMES), len(data)))

		assert container == container2

		# Check container2's internals after custom unpickling.
		np.testing.assert_array_equal(container._counts, container2._counts)
		assert container._nObjects == container2._nObjects
		assert container._dtype == container2._dtype
		assert container._objectIndex == container2._objectIndex

		# Test an int16 container.
		container3 = BulkObjectsContainer(ELEMENTS, dtype=np.int16)
		container3.countsIs(301 * np.arange(len(ELEMENTS)))

		data = cPickle.dumps(container3, cPickle.HIGHEST_PROTOCOL)
		container4 = cPickle.loads(data)
		# print("Pickled a BulkObjectsContainer of {} int16 to {} bytes".format(
		# 	len(ELEMENT_NAMES), len(data)))

		assert container4 == container3

		# Test that an unpickled container is mutable as always. It's at risk
		# because np.frombuffer(f.read()) gets an ndarray onto a readonly byte str.
		e = container4.count('Einsteinium')
		container4.countInc(1010, 'Einsteinium')
		self.assertEqual(1010 + e, container4.count('Einsteinium'))

	def test_cannot_pickle(self):
		"""Try to pickle a container whose dtype has fields or a subarray."""
		container = BulkObjectsContainer(ELEMENTS, dtype=[('U', 'f4'), ('b', 'i4')])
		with self.assertRaises(ValueError):
			cPickle.dumps(container)

		container = BulkObjectsContainer(ELEMENTS, dtype='(2,3)f8')
		with self.assertRaises(ValueError):
			cPickle.dumps(container)

		with self.assertRaises(TypeError):
			# Suppress PyCharm's type check then test that the bad arg type
			# gets caught at run time -- currently by BOC's dumps() method.
			# Why doesn't mypy catch this bad type?
			# noinspection PyTypeChecker
			bad_names = [OBJECT_NAMES, OBJECT_NAMES]  # type: Iterable[str]
			container = BulkObjectsContainer(bad_names, dtype='f8')
			cPickle.dumps(container)

	def test_write_table(self):
		"""Test writing the container to a Table."""
		self.make_test_dir()
		path = os.path.join(self.test_dir, 'BulkMolecules')

		container = BulkObjectsContainer(ELEMENTS, dtype=np.float64)
		container.countsIs(234 * np.arange(len(ELEMENTS)))  # values fit in in16

		table_writer = TableWriter(path)
		container.tableCreate(table_writer)
		container.tableAppend(table_writer)
		table_writer.close()

		# Read and check the data.
		# ASSUMES: One table row appended, readColumn() squeezes the array to
		# 1D, and various internal details about how this container writes a
		# table such as keeping unused array elements.
		table_reader = TableReader(path)
		self.assertEqual(ELEMENTS, table_reader.readAttribute('objectNames'))
		npt.assert_array_equal(container.counts(), table_reader.readColumn('counts'))


# TODO: view tests


def createContainer():
	container = BulkObjectsContainer(OBJECT_NAMES)

	container.countsIs(OBJECT_COUNTS)

	return container
