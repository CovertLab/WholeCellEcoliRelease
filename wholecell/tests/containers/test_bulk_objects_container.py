'''
test_bulk_objects_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

import wholecell.utils.bulk_objects_container

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
		self.assertEqual(
			self.container.counts().tolist(),
			OBJECT_COUNTS
			)

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			OBJECT_COUNTS[1:]
			)

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsIs(self):
		newCounts = [10, 20, 30]
		self.container.countsIs(newCounts)

		self.assertEqual(
			self.container.counts().tolist(),
			newCounts
			)

		newCounts = [15, 5]
		self.container.countsIs(newCounts, OBJECT_NAMES[1:])

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			newCounts
			)

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsInc(self):
		incCounts = [10, 20, 30]
		newCounts = (self.container.counts() + np.array(incCounts)).tolist()
		self.container.countsInc(incCounts)

		self.assertEqual(
			self.container.counts().tolist(),
			newCounts
			)

		incCounts = [15, 5]
		newCounts = (self.container.counts(OBJECT_NAMES[1:]) + np.array(incCounts)).tolist()
		self.container.countsInc(incCounts, OBJECT_NAMES[1:])

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			newCounts
			)

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsDec(self):
		decCounts = [30, 10, 5]
		newCounts = (np.array(OBJECT_COUNTS) - np.array(decCounts)).tolist()
		self.container.countsDec(decCounts)

		self.assertEqual(
			self.container.counts().tolist(),
			newCounts
			)

		decCounts = [5, 2]
		newCounts = (self.container.counts(OBJECT_NAMES[1:]) - np.array(decCounts)).tolist()
		self.container.countsDec(decCounts, OBJECT_NAMES[1:])

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			newCounts
			)

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_count(self):
		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			OBJECT_COUNTS[0]
			)

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countIs(self):
		newCount = 40
		self.container.countIs(newCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countInc(self):
		incCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countDec(self):
		decCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)

	# Internal methods

	
	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_namesToIndexes(self):
		# Test normal ordering
		names = OBJECT_NAMES
		self.assertEqual(
			self.container._namesToIndexes(names).tolist(),
			np.arange(len(OBJECT_NAMES)).tolist()
			)

		# Test reverse ordering
		names = OBJECT_NAMES[::-1]
		self.assertEqual(
			self.container._namesToIndexes(names).tolist(),
			np.arange(len(OBJECT_NAMES))[::-1].tolist()
			)

		# Test subset
		names = OBJECT_NAMES[::2]
		self.assertEqual(
			self.container._namesToIndexes(names).tolist(),
			np.arange(len(OBJECT_NAMES))[::2].tolist()
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_eq(self):
		newContainer = createContainer()

		self.assertEqual(
			self.container.counts().tolist(),
			newContainer.counts().tolist()
			)

# TODO: view tests


def createContainer():
	container = wholecell.utils.bulk_objects_container.BulkObjectsContainer(
		OBJECT_NAMES)

	container.countsIs(OBJECT_COUNTS)

	return container
