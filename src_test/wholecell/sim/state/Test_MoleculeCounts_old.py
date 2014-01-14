#!/usr/bin/env python

"""
Test MoleculeCounts.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/5/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools

import numpy
import cPickle
import os
# import wholecell.util.randStream
import wholecell.sim.state.MoleculeCounts


class Test_MoleculeCounts(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = None
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		self.moleculeCounts = self.sim.getState('MoleculeCounts')

		# Create generic process for partition
		self.genericProcess = type("", (), {})()
		self.genericProcess.meta = {"id": "genericProcess_id", "name": "genericProcess_name"}

		# Set counts for partitioning
		self.metaboliteIdx = self.moleculeCounts.getIndex(["ATP[c]","CTP[c]","GTP[c]","UTP[c]","PPI[c]","H2O[c]","H[c]"])[1]
		self.moleculeCounts.counts = numpy.zeros(self.moleculeCounts.counts.shape)
		self.moleculeCounts.counts[self.metaboliteIdx,0] = numpy.array([10.,2.,5.,7.,20.,3.,7.])
		
		# Re-write partitions
		self.moleculeCounts.partitions = []

	def tearDown(self):
		pass

	@noseAttrib.attr('partitionTest')
	def test_relativeAllocation(self):
		self.metabolitePartition1 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.metabolitePartition2 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.metabolitePartition3 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.moleculeCounts.prepartition()
		self.moleculeCounts.partition()
		
		self.assertEqual(self.moleculeCounts.partitions[0].counts.tolist(), [1., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.moleculeCounts.partitions[1].counts.tolist(), [1., 1., 2., 4., 0., 0., 0.])
		self.assertEqual(self.moleculeCounts.partitions[2].counts.tolist(), [7., 0., 3., 2., 20., 0., 7.])

	@noseAttrib.attr('partitionTest')
	def test_absoluteAllocation(self):
		self.metabolitePartition1 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.metabolitePartition2 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([5., 3., 2., 2., 0., 0., 0.]), isReqAbs = True)
		self.metabolitePartition3 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([20., 1., 3., 1., 2., 0., 2.]))

		self.moleculeCounts.prepartition()
		self.moleculeCounts.partition()
		
		self.assertEqual(self.moleculeCounts.partitions[0].counts.tolist(), [0., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.moleculeCounts.partitions[1].counts.tolist(), [5., 2., 2., 2., 0., 0., 0.])
		self.assertEqual(self.moleculeCounts.partitions[2].counts.tolist(), [4., 0., 3., 5., 20., 0., 7.])


	@noseAttrib.attr('partitionTest')
	def test_absoluteAllocation_withConflict(self):
		self.metabolitePartition1 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([3., 0., 0., 0., 0., 0., 0.]), isReqAbs = False)
		self.metabolitePartition2 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([5., 3., 2., 2., 0., 0., 0.]), isReqAbs = True)
		self.metabolitePartition3 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], lambda: numpy.array([20., 1., 3., 1., 2., 0., 2.]), isReqAbs = True)

		self.moleculeCounts.prepartition()
		self.moleculeCounts.partition()
		
		self.assertEqual(self.moleculeCounts.partitions[0].counts.tolist(), [0., 0., 0., 0., 0., 0., 0.])
		self.assertEqual(self.moleculeCounts.partitions[1].counts.tolist(), [2., 1., 2., 2., 0., 0., 0.])
		self.assertEqual(self.moleculeCounts.partitions[2].counts.tolist(), [8., 0., 3., 1., 2., 0., 2.])
