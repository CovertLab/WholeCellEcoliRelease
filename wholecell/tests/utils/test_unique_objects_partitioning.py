'''
test_unique_objects_partitioning.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/4/2014
'''

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

import wholecell.utils.unique_objects_container
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
		self.container = wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			TEST_KB)
		
		self.container.moleculesNew(
			'A',
			60,
			attribute = True,
			)

		self.container.moleculesNew(
			'B',
			20,
			attribute = True,
			)

		self.container.moleculesNew(
			'B',
			20,
			attribute = False,
			)
		
		self.arrayIndex_A = self.container._nameToArrayIndex['A']
		self.arrayIndex_B = self.container._nameToArrayIndex['B']

		self.randStream = wholecell.utils.rand_stream.RandStream()


	def tearDown(self):
		pass

	
	@noseAttrib.attr('smalltest', 'uniqueObjects', 'partitioning')
	def test_partitioning(self):
		request11_MoleculesLocal = self.container._queryMolecules(self.arrayIndex_A,
			attribute = ('==', True))

		request12_MoleculesLocal = self.container._queryMolecules(self.arrayIndex_B,
			attribute = ('==', True))

		request21_MoleculesLocal = self.container._queryMolecules(self.arrayIndex_A,
			attribute = ('==', True))

		request22_MoleculesLocal = self.container._queryMolecules(self.arrayIndex_B)

		globalRefSize = self.container._arrays[self.container._globalRefIndex].size

		globalIndexes_A = self.container._arrays[self.arrayIndex_A]['_globalIndex']
		globalIndexes_B = self.container._arrays[self.arrayIndex_B]['_globalIndex']

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

		wholecell.utils.unique_objects_container._partition(objectRequestsArray,
			requestNumberVector, requestProcessArray, self.randStream)


	
