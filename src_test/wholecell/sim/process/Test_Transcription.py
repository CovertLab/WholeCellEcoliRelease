"""
Test Transcription.py
Tests Transcription process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/16/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

from mpi4py import MPI

import numpy
import scipy.stats
import cPickle
import os
import matplotlib
matplotlib.use("agg")
from wholecell.util.Constants import Constants

comm = MPI.COMM_WORLD

import copy_reg
import types

# We need the following two methods so that we can pickle instancemethods
# Also need the copy_reg.pickle() line in generateTestFixtures()
# Borrowed from: http://mail.python.org/pipermail/python-list/2006-October/367078.html
# (Thanks Steven Bethard!)

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

class Test_Transcription(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = None

		if comm.rank == 0:
			self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		
		self.sim = comm.bcast(self.sim, root = 0)	
		print "%s" % (self.sim.getState("Mass").meta["id"])

	def tearDown(self):
		pass


	# Tests

	def test_production(self):
		sim = self.sim
		tc = sim.getProcess("Transcription")
		tc.rna.mws[tc.rna.mws < 0 ] = 0
		
		ntpCounts = 1e6
		initEnzCnts = 2000.
		initRnaCnts = 0.
		T_d = 3600.
		lengthSec = 100

		nSeeds = 10
		nProc = comm.size
		sendcounts = (nSeeds / nProc) * numpy.ones(nProc, dtype = int) + (numpy.arange(nProc, dtype = int) < nSeeds % nProc)
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(sendcounts)[:-1]])
		mySeeds = numpy.zeros(sendcounts[comm.rank])
		allSeeds = numpy.arange(nSeeds, dtype = float)
		comm.Scatterv([allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE], mySeeds)
		print "Rank %d mySeeds: %s" % (comm.rank, str(mySeeds))

		myNtpUsage = numpy.zeros((tc.metabolite.idx["ntps"].size, lengthSec, mySeeds.size))
		myRnaProduction = numpy.zeros((tc.rna.counts.size, lengthSec, mySeeds.size))

		for iterSeed in xrange(mySeeds.size):
			seed = mySeeds.astype("int")[iterSeed]
			print "%d" % seed
			tc.randStream.seed = seed
			tc.rna.counts = initRnaCnts * numpy.ones(tc.rna.counts.shape)

			ntpUsage = numpy.zeros((tc.metabolite.idx["ntps"].size, lengthSec))
			rnaProduction = numpy.zeros((tc.rna.counts.size, lengthSec))

			for t in xrange(lengthSec):
				tc.metabolite.counts[tc.metabolite.idx["ntps"]] = ntpCounts * numpy.ones(tc.metabolite.idx["ntps"].shape)
				tc.enzyme.counts = numpy.round(initEnzCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tc.enzyme.counts.shape)
				tc.evolveState()
				ntpUsage[:, t] = tc.metabolite.parentState.tcNtpUsage
				rnaProduction[:, t] = tc.rna.counts

			myNtpUsage[:, :, iterSeed] = ntpUsage
			myRnaProduction[:, :, iterSeed] = rnaProduction

			fNtpUsage = ntpUsage / numpy.sum(ntpUsage, axis = 0)
			print "%d, %s" % (comm.rank, str(numpy.mean(fNtpUsage, axis = 1)))
			rnaMassProduction = numpy.diff(numpy.dot(tc.rna.mws, rnaProduction / Constants.nAvogadro))

			# Assert exponential usage of metabolites
			self.assertTrue(numpy.allclose(1., numpy.mean(numpy.mean(ntpUsage[:, -10:], axis = 1) / numpy.mean(ntpUsage[:, :10], axis = 1) / numpy.exp(numpy.log(2) / T_d * lengthSec)), rtol = 0, atol = 6e-2))
			self.assertTrue(numpy.allclose(1., numpy.mean(ntpUsage[:, -10:], axis = 1) / numpy.mean(ntpUsage[:, :10], axis = 1) / numpy.exp(numpy.log(2) / T_d * lengthSec), rtol = 0, atol = 6e-2))

			# Assert exponential increase in RNA mass
			self.assertTrue(numpy.allclose(1., numpy.mean(rnaMassProduction[-10:]) / numpy.mean(rnaMassProduction[:10]) / numpy.exp(numpy.log(2) / T_d * lengthSec), rtol = 0, atol = 6e-2))

			nProd = rnaProduction[:, -1]
			N = numpy.sum(nProd)
			E = N * tc.rnaSynthProb
			V = N * tc.rnaSynthProb * (1 - tc.rnaSynthProb)
			S = numpy.sqrt(V)

			# Assert exponential increase in rRNA 23S counts
			rrlIdxs = tc.rna.getIndex(["RRLA-RRNA:nascent", "RRLB-RRNA:nascent", "RRLC-RRNA:nascent", "RRLD-RRNA:nascent", "RRLE-RRNA:nascent", "RRLG-RRNA:nascent", "RRLH-RRNA:nascent"])[0]
			self.assertTrue(numpy.allclose(1., numpy.mean(nProd[rrlIdxs] / E[rrlIdxs]), rtol = 0, atol = 1e-1))

			# Loosely assert exponential increase in all RNAs (this is noisy)
			numpyHandleError = numpy.geterr()
			numpy.seterr(divide = "ignore", invalid = "ignore")
			self.assertTrue(scipy.stats.nanmean(nProd / E) > 0.8)
			numpy.seterr(**numpyHandleError)

		if comm.rank == 0:
			allNtpUsage = numpy.zeros(tc.metabolite.idx["ntps"].size * lengthSec * allSeeds.size)
			allRnaProduction = numpy.zeros(tc.rna.counts.size * lengthSec * allSeeds.size)
		else:
			allNtpUsage = None
			allRnaProduction = None

		recvcounts = sendcounts * ntpUsage.size
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
		comm.Gatherv(numpy.concatenate(numpy.split(myNtpUsage, mySeeds.size, axis = 2)), [allNtpUsage, recvcounts, displacements, MPI.DOUBLE])
		recvcounts = sendcounts * rnaProduction.size
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
		comm.Gatherv(numpy.concatenate(numpy.split(myRnaProduction, mySeeds.size, axis = 2)), [allRnaProduction, recvcounts, displacements, MPI.DOUBLE])

		if comm.rank == 0:
			allNtpUsage = numpy.dstack([x.reshape(ntpUsage.shape) for x in numpy.split(allNtpUsage.reshape(-1), allSeeds.size)])
			allRnaProduction = numpy.dstack([x.reshape(rnaProduction.shape) for x in numpy.split(allRnaProduction.reshape(-1), allSeeds.size)])
			print "%f == %f" % (allRnaProduction[10, 10, 0], myRnaProduction[10, 10, 0])
			print "%f == %f" % (allRnaProduction[10, 50, 3], myRnaProduction[10, 50, 3])
			print "%f == %f" % (allRnaProduction[-1, 50, 2], myRnaProduction[-1, 50, 2])
			import ipdb
			ipdb.set_trace()
			numpy.all(allRnaProduction[:, :, 0:5] == myRnaProduction)