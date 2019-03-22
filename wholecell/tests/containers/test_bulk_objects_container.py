'''
test_bulk_objects_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import absolute_import
from __future__ import division

import unittest

import numpy as np
import numpy.testing as npt
import nose.plugins.attrib as noseAttrib

from wholecell.containers.bulk_objects_container import BulkObjectsContainer

OBJECT_NAMES = ('ATP', 'glucose', 'glycine')
OBJECT_COUNTS = [100, 20, 10]


class Test_BulkObjectsContainer(unittest.TestCase):
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

	# Interface methods

	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_counts(self):
		npt.assert_equal(self.container.counts(), OBJECT_COUNTS)

		npt.assert_equal(
			self.container.counts(OBJECT_NAMES[1:]),
			OBJECT_COUNTS[1:]
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsIs(self):
		newCounts = [10, 20.0, 30]
		self.container.countsIs(newCounts)

		npt.assert_equal(self.container.counts(), newCounts)

		newCounts = [15, 5.0]
		self.container.countsIs(newCounts, OBJECT_NAMES[1:])

		npt.assert_equal(self.container.counts(OBJECT_NAMES[1:]), newCounts)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsInc(self):
		incCounts = [10, 20.0, 30]
		newCounts = self.container.counts() + incCounts
		self.container.countsInc(incCounts)

		npt.assert_equal(self.container.counts(), newCounts)

		incCounts = [15, 5.0]
		newCounts = self.container.counts(OBJECT_NAMES[1:]) + incCounts
		self.container.countsInc(incCounts, OBJECT_NAMES[1:])

		npt.assert_equal(self.container.counts(OBJECT_NAMES[1:]), newCounts)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsDec(self):
		decCounts = [30, 10.0, 5]
		newCounts = np.array(OBJECT_COUNTS) - np.array(decCounts)
		self.container.countsDec(decCounts)

		npt.assert_equal(self.container.counts(), newCounts)

		decCounts = [5.0, 2]
		newCounts = self.container.counts(OBJECT_NAMES[1:]) - np.array(decCounts)
		self.container.countsDec(decCounts, OBJECT_NAMES[1:])

		npt.assert_equal(self.container.counts(OBJECT_NAMES[1:]), newCounts)


	@noseAttrib.attr('smalltest', 'bulkObjects')
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


	@noseAttrib.attr('smalltest', 'bulkObjects')
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


	@noseAttrib.attr('smalltest', 'bulkObjects')
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


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countView(self):
		countView = self.container.countView(OBJECT_NAMES[2])
		incCount = 333.0
		newCount = OBJECT_COUNTS[2] + incCount
		countView.countInc(incCount)

		npt.assert_equal(countView.count(), newCount)

		newCount -= incCount
		countView.countDec(incCount)
		npt.assert_equal(countView.count(), newCount)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_count(self):
		npt.assert_equal(
			self.container.count(OBJECT_NAMES[0]),
			OBJECT_COUNTS[0]
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countIs(self):
		newCount = 40
		self.container.countIs(newCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countInc(self):
		incCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)

		incCount = 40.0
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countDec(self):
		decCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)

		decCount = 40.0
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		npt.assert_equal(self.container.count(OBJECT_NAMES[0]), newCount)


	# Internal methods


	@noseAttrib.attr('smalltest', 'bulkObjects')
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


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_eq(self):
		newContainer = createContainer()

		npt.assert_equal(self.container.counts(), newContainer.counts())

# TODO: view tests


def createContainer():
	container = BulkObjectsContainer(OBJECT_NAMES)

	container.countsIs(OBJECT_COUNTS)

	return container
