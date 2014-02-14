"""
Test RnaDegradation.py
Tests RNA degradation process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/26/2013
"""

# import unittest
# import warnings
# import nose.plugins.attrib as noseAttrib

# import numpy
# import scipy.stats
# import cPickle
# import os
# #import matplotlib
# #matplotlib.use("agg")
# from wholecell.util.Constants import Constants

# from mpi4py import MPI

# comm = MPI.COMM_WORLD

# import copy_reg
# import types

# # We need the following two methods so that we can pickle instancemethods
# # Also need the copy_reg.pickle() line in generateTestFixtures()
# # Borrowed from: http://mail.python.org/pipermail/python-list/2006-October/367078.html
# # (Thanks Steven Bethard!)

# def _pickle_method(method):
#     func_name = method.im_func.__name__
#     obj = method.im_self
#     cls = method.im_class
#     return _unpickle_method, (func_name, obj, cls)

# def _unpickle_method(func_name, obj, cls):
#     for cls in cls.mro():
#         try:
#             func = cls.__dict__[func_name]
#         except KeyError:
#             pass
#         else:
#             break
#     return func.__get__(obj, cls)

# copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

# class Test_RnaDegradation(unittest.TestCase):

# 	@classmethod
# 	def setUpClass(cls):
# 		pass

# 	@classmethod
# 	def tearDownClass(cls):
# 		pass

# 	def setUp(self):
# 		self.sim = None

# 		if comm.rank == 0:
# 			self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		
# 		self.sim = comm.bcast(self.sim, root = 0)	
# 		print "%s" % (self.sim.states["Mass"].meta["id"])

# 	def tearDown(self):
# 		pass


# 	# Tests

# 	# TODO: Parallelize this test
# 	@noseAttrib.attr('largetest')
# 	def test_degradation(self):
# 		sim = self.sim
# 		rd = sim.processes["RnaDegradation"]
# 		mc = sim.states["MoleculeCounts"]
# 		rd.rna.mws[rd.rna.mws < 0 ] = 0
		
# 		h2oCounts = 1e6
# 		initEnzCnts = 1000.
# 		initRnaCnts = 10000. * numpy.ones(rd.rna.counts.shape)
# 		T_d = 3600.
# 		lengthSec = 100

# 		for iterSeed in xrange(100):
# 			rd.randStream.seed = iterSeed

# 			rd.rna.fullCounts = initRnaCnts.copy()
# 			nDeg = numpy.zeros((rd.rna.counts.size, lengthSec))

# 			for t in xrange(lengthSec):
# 				rd.rna.counts = rd.calcReqRna()
# 				begRnaCounts = rd.rna.counts.copy()
# 				rd.rna.fullCounts -= rd.rna.counts
# 				rd.metabolite.counts = numpy.zeros(rd.metabolite.counts.shape)
# 				rd.metabolite.counts[rd.metabolite.idx["h2o"]] = h2oCounts
# 				rd.enzyme.counts[:] = initEnzCnts
# 				rd.evolveState()
# 				rd.rna.fullCounts += rd.rna.counts

# 				nDeg[:, t] += (begRnaCounts - rd.rna.counts)

# 			# Assert that we accounted for degradation (nDeg) properly
# 			self.assertTrue(numpy.all(rd.rna.fullCounts + numpy.sum(nDeg, axis = 1) == initRnaCnts))

# 			# Assert that we have approximately exponential decay
# 			self.assertTrue(numpy.max(numpy.abs((rd.rna.fullCounts - numpy.exp(-1 * rd.rnaDegRates * lengthSec) * initRnaCnts)) / (numpy.exp(-1 * rd.rnaDegRates * lengthSec) * initRnaCnts)) < 5e-2)
		
# 			# TODO: Assert mass balance (look at metabolite counts and stuff)
# 			# Makes sense to do when Metabolism is modeled better

