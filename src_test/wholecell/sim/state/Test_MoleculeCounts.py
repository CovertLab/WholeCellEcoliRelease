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
		self.metaboliteIdx = self.moleculeCounts.getIndex(["ATP[c]","CTP[c]","GTP[c]","UTP[c]","PPI[c]","H2O[c]","H[c]"])[0]
		self.moleculeCounts.counts = numpy.zeros(self.moleculeCounts.counts.shape)
		self.moleculeCounts.counts[self.metaboliteIdx,0] = numpy.array([10.,2.,5.,7.,20.,3.,7.])
		
		# Re-write partitions
		self.moleculeCounts.partitions = []
		self.metabolitePartition1 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], None)
		self.metabolitePartition1.reqFunc = make_testRequestFunction(self.metabolitePartition1,numpy.array([3., 0., 0., 0., 0., 0., 0.]))
		self.metabolitePartition2 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], None)
		self.metabolitePartition2.reqFunc = make_testRequestFunction(self.metabolitePartition2,numpy.array([5., 3., 2., 2., 0., 0., 0.]))
		self.metabolitePartition3 = self.moleculeCounts.addPartition(self.genericProcess, [
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"PPI[c]", "H2O[c]", "H[c]",
			], None)
		self.metabolitePartition3.reqFunc = make_testRequestFunction(self.metabolitePartition3,numpy.array([20., 1., 3., 1., 2., 0., 2.]))


	def tearDown(self):
		pass

	@noseAttrib.attr('partitionTest')
	def test_enzymeAllocation(self):
		self.moleculeCounts.prepartition()
		self.moleculeCounts.partition()
		import ipdb; ipdb.set_trace()

		#self.assertEqual(enzymeAlloc_transcription.counts.tolist(), numpy.array([  958.,  1000.,  1014.,   980.]).tolist())

		# "EG10893-MONOMER", "RPOB-MONOMER", "RPOC-MONOMER", "RPOD-MONOMER"

def make_testRequestFunction(partition, request):
	def testRequestFunction():
		return numpy.array(request)
	return testRequestFunction