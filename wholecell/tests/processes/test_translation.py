"""
Test Translation.py
Tests Translation process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/16/2013
"""

# import unittest
# import warnings
# import nose.plugins.attrib as noseAttrib

# import numpy
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

# class Test_Translation(unittest.TestCase):

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
# 	@noseAttrib.attr('largetest')
# 	def test_production(self):

# 		sim = self.sim
# 		tl = sim.processes["Translation"]
# 		tl.protein.mws[tl.protein.mws < 0 ] = 0
		
# 		aaCounts = 1e6
# 		initEnzCnts = 1000.
# 		initProtCnts = 200.
# 		initRnaCnts = 100.
# 		T_d = 3600.
# 		lengthSec = 1000
# 		notSecIdxs = tl.metabolite.idx["aasNotSec"]

# 		nSeeds = 10
# 		nProc = comm.size
# 		sendcounts = (nSeeds / nProc) * numpy.ones(nProc, dtype = int) + (numpy.arange(nProc, dtype = int) < nSeeds % nProc)
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(sendcounts)[:-1]])
# 		mySeeds = numpy.zeros(sendcounts[comm.rank])
# 		allSeeds = numpy.arange(nSeeds, dtype = float)
# 		comm.Scatterv([allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE], mySeeds)
# 		print "Rank %d mySeeds: %s" % (comm.rank, str(mySeeds))

# 		myAaUsage = numpy.zeros((tl.metabolite.idx["aas"].size, lengthSec, mySeeds.size))
# 		myMonProduction = numpy.zeros((tl.protein.counts.size, lengthSec, mySeeds.size))

# 		for iterSeed in xrange(mySeeds.size):
# 			seed = mySeeds.astype("int")[iterSeed]
# 			tl.randStream.seed = seed
# 			tl.protein.counts = initProtCnts * numpy.ones(tl.protein.counts.shape)

# 			aaUsage = numpy.zeros((tl.metabolite.idx["aas"].size, lengthSec))
# 			monProduction = numpy.zeros((tl.protein.counts.size, lengthSec))

# 			for t in xrange(lengthSec):
# 				tl.metabolite.counts[tl.metabolite.idx["aas"]] = aaCounts * numpy.ones(tl.metabolite.idx["aas"].shape)
# 				tl.enzyme.counts = numpy.round(initEnzCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tl.enzyme.counts.shape)
# 				tl.mrna.counts = numpy.round(initRnaCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tl.mrna.counts.shape)
# 				tl.evolveState()
# 				tl.aasUsed[15] += 1e-9
# 				aaUsage[:, t] = tl.aasUsed
# 				monProduction[:, t] = tl.protein.counts

# 			myAaUsage[:, :, iterSeed] = aaUsage
# 			myMonProduction[:, :, iterSeed] = monProduction

# 		## MPI communications ##
# 		# Get ready to gather values
# 		if comm.rank == 0:
# 			allAaUsage = numpy.zeros(tl.metabolite.idx["aas"].size * lengthSec * allSeeds.size)
# 			allMonProduction = numpy.zeros(tl.protein.counts.size * lengthSec * allSeeds.size)
# 		else:
# 			allAaUsage = None
# 			allMonProduction = None

# 		# Gather AA usage values
# 		recvcounts = sendcounts * aaUsage.size
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myAaUsage, mySeeds.size, axis = 2)), [allAaUsage, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather monomer production values
# 		recvcounts = sendcounts * monProduction.size
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myMonProduction, mySeeds.size, axis = 2)), [allMonProduction, recvcounts, displacements, MPI.DOUBLE])

# 		# Assemble gathered values
# 		if comm.rank == 0:
# 			allAaUsage = numpy.dstack([x.reshape(aaUsage.shape) for x in numpy.split(allAaUsage.reshape(-1), allSeeds.size)])
# 			allMonProduction = numpy.dstack([x.reshape(myMonProduction.shape) for x in numpy.split(allMonProduction.reshape(-1), allSeeds.size)])
			
# 			# Check for proper assembly of the gathered data
# 			# (Check that it's not incorrect, not necessarily guaranteed to be correct)
# 			self.assertTrue(numpy.all(allAaUsage[:, :, :mySeeds.size] == myAaUsage))
# 			self.assertTrue(numpy.all(allMonProduction[:, :, :mySeeds.size] == myMonProduction))

# 			print "===== Amino acid usage ====="
# 			for iDepth in xrange(allAaUsage.shape[2]):
# 				thisAaUsage = numpy.squeeze(allAaUsage[:, :, iDepth])
# 				print "[%d] Mean of each AA relative to expectation: %0.3f" % (iDepth, numpy.mean(numpy.mean(thisAaUsage[notSecIdxs, -10:], axis = 1) / numpy.mean(thisAaUsage[notSecIdxs, :10], axis = 1) / numpy.exp(numpy.log(2) / T_d * lengthSec)))
# 				self.assertTrue(numpy.allclose(1., numpy.mean(numpy.mean(thisAaUsage[notSecIdxs, -10:], axis = 1) / numpy.mean(thisAaUsage[notSecIdxs, :10], axis = 1) / numpy.exp(numpy.log(2) / T_d * lengthSec)), rtol = 0, atol = 1e-1))
# 				self.assertTrue(numpy.allclose(1., numpy.mean(thisAaUsage[notSecIdxs, -10:], axis = 1) / numpy.mean(thisAaUsage[notSecIdxs, :10], axis = 1) / numpy.exp(numpy.log(2) / T_d * lengthSec), rtol = 0, atol = 12e-2))

# 			print "===== Monomer mass production ====="
# 			for iDepth in xrange(allMonProduction.shape[2]):
# 				thisMonProduction = numpy.squeeze(allMonProduction[:, :, iDepth])
# 				thisMonMassProduction = numpy.diff(numpy.dot(tl.protein.mws, thisMonProduction / Constants.nAvogadro))
# 				self.assertTrue(numpy.mean(thisMonMassProduction[-10:]) / numpy.mean(thisMonMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec) > 0.98)
# 				self.assertTrue(numpy.mean(thisMonMassProduction[-10:]) / numpy.mean(thisMonMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec) < 1.08)
# 				print "[%d] Mass production relative to expectation: %0.3f" % (iDepth, numpy.mean(thisMonMassProduction[-10:]) / numpy.mean(thisMonMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec))

# 			# TODO: More assertions
# 			# - Exponential increase in production of all monomers
			

# 	@noseAttrib.attr("monTotalProduction")
# 	@noseAttrib.attr('largetest')
# 	def test_total_production(self):
# 		sim = self.sim
# 		tl = sim.processes["Translation"]
# 		tl.protein.mws[tl.protein.mws < 0 ] = 0
# 		mc = sim.states["MoleculeCounts"]
		
# 		aaCounts = 1e6
# 		initEnzCnts = 1100.
# 		initProtCnts = 0.
# 		T_d = 3600.
# 		lengthSec = int(T_d)

# 		nSeeds = 10
# 		nProc = comm.size
# 		sendcounts = (nSeeds / nProc) * numpy.ones(nProc, dtype = int) + (numpy.arange(nProc, dtype = int) < nSeeds % nProc)
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(sendcounts)[:-1]])
# 		mySeeds = numpy.zeros(sendcounts[comm.rank])
# 		allSeeds = numpy.arange(nSeeds, dtype = float)
# 		comm.Scatterv([allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE], mySeeds)
# 		print "Rank %d mySeeds: %s" % (comm.rank, str(mySeeds))

# 		myAaUsage = numpy.zeros((tl.metabolite.idx["aas"].size, lengthSec, mySeeds.size))
# 		myMonProduction = numpy.zeros((tl.protein.counts.size, lengthSec, mySeeds.size))

# 		for iterSeed in xrange(mySeeds.size):
# 			seed = mySeeds.astype("int")[iterSeed]
# 			tl.randStream.seed = seed
# 			tl.protein.counts = initProtCnts * numpy.ones(tl.protein.counts.shape)

# 			mc.randStream.seed = seed
# 			mc.calcInitialConditions()
# 			initRnaCnts = mc.counts[mc.idx["matureMrna"], mc.cIdx["c"]]

# 			aaUsage = numpy.zeros((tl.metabolite.idx["aas"].size, lengthSec))
# 			monProduction = numpy.zeros((tl.protein.counts.size, lengthSec))

# 			for t in xrange(lengthSec):
# 				tl.metabolite.counts[tl.metabolite.idx["aas"]] = aaCounts * numpy.ones(tl.metabolite.idx["aas"].shape)
# 				tl.enzyme.counts = numpy.round(initEnzCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tl.enzyme.counts.shape)
# 				tl.mrna.counts = numpy.round(initRnaCnts * numpy.exp(numpy.log(2) / T_d * t)) #* numpy.ones(tl.mrna.counts.shape)
# 				tl.evolveState()
# 				tl.aasUsed[15] += 1e-9
# 				aaUsage[:, t] = tl.aasUsed
# 				monProduction[:, t] = tl.protein.counts

# 			myAaUsage[:, :, iterSeed] = aaUsage
# 			myMonProduction[:, :, iterSeed] = monProduction

# 		## MPI communications ##
# 		# Get ready to gather values
# 		if comm.rank == 0:
# 			allAaUsage = numpy.zeros(tl.metabolite.idx["aas"].size * lengthSec * allSeeds.size)
# 			allMonProduction = numpy.zeros(tl.protein.counts.size * lengthSec * allSeeds.size)
# 		else:
# 			allAaUsage = None
# 			allMonProduction = None

# 		# Gather AA usage values
# 		recvcounts = sendcounts * aaUsage.size
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myAaUsage, mySeeds.size, axis = 2)), [allAaUsage, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather monomer production values
# 		recvcounts = sendcounts * monProduction.size
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myMonProduction, mySeeds.size, axis = 2)), [allMonProduction, recvcounts, displacements, MPI.DOUBLE])

# 		# Assemble gathered values
# 		if comm.rank == 0:
# 			allAaUsage = numpy.dstack([x.reshape(aaUsage.shape) for x in numpy.split(allAaUsage.reshape(-1), allSeeds.size)])
# 			allMonProduction = numpy.dstack([x.reshape(myMonProduction.shape) for x in numpy.split(allMonProduction.reshape(-1), allSeeds.size)])
			
# 			# Check for proper assembly of the gathered data
# 			# (Check that it's not incorrect, not necessarily guaranteed to be correct)
# 			self.assertTrue(numpy.all(allAaUsage[:, :, :mySeeds.size] == myAaUsage))
# 			self.assertTrue(numpy.all(allMonProduction[:, :, :mySeeds.size] == myMonProduction))

# 			# Assert that cumulative production is a cell's worth of protein
# 			for iDepth in xrange(allMonProduction.shape[2]):
# 				thisMonProduction = numpy.squeeze(allMonProduction[:, :, iDepth])
# 				massProd = numpy.dot(tl.protein.mws, thisMonProduction[:, -1] / Constants.nAvogadro)
# 				self.assertTrue(numpy.allclose(1., massProd / (2.8e-13 * 0.563 / 1.36), rtol = 0, atol = 2e-2))
