'''
test_unique_objects_partitioning.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/4/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.containers.unique_objects_container import UniqueObjectsContainer, _partition
import wholecell.utils.rand_stream

TEST_KB = {
	'A':{
		'attribute':'bool',
		},
	'B':{
		'attribute':'bool',
		}
	}


class Test_UniqueMoleculesContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = UniqueObjectsContainer(TEST_KB)
		
		self.container.objectsNew(
			'A',
			60,
			attribute = True,
			)

		self.container.objectsNew(
			'B',
			20,
			attribute = True,
			)

		self.container.objectsNew(
			'B',
			20,
			attribute = False,
			)
		
		self.arrayIndex_A = self.container._nameToCollectionIndex['A']
		self.arrayIndex_B = self.container._nameToCollectionIndex['B']

		self.randStream = wholecell.utils.rand_stream.RandStream()


	def tearDown(self):
		pass

	
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'partitioning')
	def test_partitioning(self):
		# Set up the partition function call
		
		# TODO: write and use the interface for this code
		request11_MoleculesLocal = self.container._queryObjects(self.arrayIndex_A,
			attribute = ('==', True))

		request12_MoleculesLocal = self.container._queryObjects(self.arrayIndex_B,
			attribute = ('==', True))

		request21_MoleculesLocal = self.container._queryObjects(self.arrayIndex_A,
			attribute = ('==', True))

		request22_MoleculesLocal = self.container._queryObjects(self.arrayIndex_B)

		globalRefSize = self.container._collections[self.container._globalRefIndex].size

		globalIndexes_A = self.container._collections[self.arrayIndex_A]['_globalIndex']
		globalIndexes_B = self.container._collections[self.arrayIndex_B]['_globalIndex']

		objectRequestsArray = np.zeros((globalRefSize, 4), np.bool)

		objectRequestsArray[globalIndexes_A, 0] = request11_MoleculesLocal
		objectRequestsArray[globalIndexes_B, 1] = request12_MoleculesLocal
		objectRequestsArray[globalIndexes_A, 2] = request21_MoleculesLocal
		objectRequestsArray[globalIndexes_B, 3] = request22_MoleculesLocal

		requestNumberVector = np.array([50, 20, 30, 20])

		requestProcessArray = np.array([
			[True, False],
			[True, False],
			[False, True],
			[False, True],
			])

		# Partition the molecules
		
		partitionedMolecules = _partition(objectRequestsArray, 
			requestNumberVector, requestProcessArray, self.randStream)

		# Assert that each molecule is partitioned to one state

		self.assertTrue((partitionedMolecules.sum(1) <= 1).all())

		# Assert that unrequested molecules aren't partitioned

		requestedByProcess1 = objectRequestsArray[:, 0] | objectRequestsArray[:, 1]

		self.assertFalse((~requestedByProcess1 & partitionedMolecules[:, 0]).any())

		# Determine the number partitioned

		partitioned11 = (objectRequestsArray[:, 0] & partitionedMolecules[:, 0]).sum()
		partitioned12 = (objectRequestsArray[:, 1] & partitionedMolecules[:, 0]).sum()
		partitioned21 = (objectRequestsArray[:, 2] & partitionedMolecules[:, 1]).sum()
		partitioned22 = (objectRequestsArray[:, 3] & partitionedMolecules[:, 1]).sum()

		# Assert that the partitioned amounts are in the correct ratios

		ratio_request1 = requestNumberVector[0]/requestNumberVector[1]
		ratio_request2 = requestNumberVector[2]/requestNumberVector[3]

		ratio_partitioned1 = partitioned11/partitioned12
		ratio_partitioned2 = partitioned21/partitioned22

		self.assertTrue(
			np.allclose(
				[ratio_request1, ratio_request2],
				[ratio_partitioned1, ratio_partitioned2]
				)
			)

	
	# def test_partitioning_worst_case_scenario(self):
	# 	nMoleculeTypes = 100
	# 	nMoleculesEach = 1000
	# 	nRequests = 100
	# 	nProcesses = 20

	# 	maxRequestCount = 3*nMoleculesEach/nRequests

	# 	objectRequestsArray = np.tile(
	# 		np.random.randint(2, size = (nMoleculeTypes, nRequests)).astype(np.bool),
	# 		(nMoleculesEach, 1)
	# 		)

	# 	requestNumberVector = np.random.randint(maxRequestCount, size = nRequests)

	# 	requestProcessArray = (np.random.randint(nProcesses, size = nRequests) == np.tile(np.arange(nProcesses), (nRequests, 1)).T).transpose()
	# 	import time; t = time.time()
	# 	solution = wholecell.utils.unique_objects_container._partition(
	# 		objectRequestsArray,requestNumberVector, requestProcessArray, self.randStream)
	# 	print time.time() - t
