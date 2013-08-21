"""
Test Translation.py
Tests Translation process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/16/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

import numpy
import cPickle
import os
import matplotlib
matplotlib.use("agg")
from wholecell.util.Constants import Constants

class Test_Translation(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))

	def tearDown(self):
		pass


	# Tests

	def test_production(self):

		sim = self.sim
		tl = sim.getProcess("Translation")
		tl.protein.mws[tl.protein.mws < 0 ] = 0
		
		aaCounts = 1e6
		initEnzCnts = 1000.
		initProtCnts = 200.
		initRnaCnts = 100.
		T_d = 3600.
		lengthSec = 1000

		for seed in xrange(10):
			tl.randStream.seed = seed
			tl.protein.counts = initProtCnts * numpy.ones(tl.protein.counts.shape)

			aaUsage = numpy.zeros((tl.metabolite.idx["aas"].size, lengthSec))
			monProduction = numpy.zeros((tl.protein.counts.size, lengthSec))

			for t in xrange(lengthSec):
				tl.metabolite.counts[tl.metabolite.idx["aas"]] = aaCounts * numpy.ones(tl.metabolite.idx["aas"].shape)
				tl.enzyme.counts = numpy.round(initEnzCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tl.enzyme.counts.shape)
				tl.mrna.counts = numpy.round(initRnaCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tl.mrna.counts.shape)
				# tl.enzyme.counts = initEnzCnts * numpy.ones(tl.enzyme.counts.shape)
				# tl.mrna.counts = initRnaCnts * numpy.ones(tl.mrna.counts.shape)
				tl.evolveState()
				aaUsage[:, t] = tl.aasUsed
				monProduction[:, t] = tl.protein.counts

			fAaUsage = aaUsage / numpy.sum(aaUsage, axis = 0)
			print "%s" % (str(numpy.mean(fAaUsage, axis = 1)))
			monMassProduction = numpy.diff(numpy.dot(tl.protein.mws, monProduction / Constants.nAvogadro))

			print "%0.3f" % (numpy.mean(numpy.mean(aaUsage[:, -10:], axis = 1) / numpy.mean(aaUsage[:, :10], axis = 1) / numpy.exp(numpy.log(2) / T_d * lengthSec)))

			print "%0.3f" % (numpy.mean(monMassProduction[-10:]) / numpy.mean(monMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec))
			self.assertTrue(numpy.mean(monMassProduction[-10:]) / numpy.mean(monMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec) > 0.95)
			self.assertTrue(numpy.mean(monMassProduction[-10:]) / numpy.mean(monMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec) < 1.10)
			# import ipdb
			# ipdb.set_trace()