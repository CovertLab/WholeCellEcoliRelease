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

import numpy
import scipy.stats
import cPickle
import os
import matplotlib
matplotlib.use("agg")
from wholecell.util.Constants import Constants

from mpi4py import MPI

comm = MPI.COMM_WORLD

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
	@noseAttrib.attr('rnaProduction')
	def test_production(self):
		sim = self.sim
		tc = sim.getProcess("Transcription")
		rnaMWs = tc.rnaPartition._state._massSingle.flat[tc.rnaPartition.mapping]
		rnaMWs[rnaMWs < 0] = 0
		
		ntpCounts = 1e6
		initEnzCnts = 2000.
		initRnaCnts = 0.
		T_d = 3600.
		lengthSec = 500

		nSeeds = 100
		nProc = comm.size
		sendcounts = (nSeeds / nProc) * numpy.ones(nProc, dtype = int) + (numpy.arange(nProc, dtype = int) < nSeeds % nProc)
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(sendcounts)[:-1]])
		mySeeds = numpy.zeros(sendcounts[comm.rank])
		allSeeds = numpy.arange(nSeeds, dtype = float)
		comm.Scatterv([allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE], mySeeds)
		print "Rank %d mySeeds: %s" % (comm.rank, str(mySeeds))

		myNtpUsage = numpy.zeros((tc.metabolitePartition.ntps.countsBulk().size, lengthSec, mySeeds.size))
		myRnaProduction = numpy.zeros((tc.rnaPartition.countsBulk().size, lengthSec, mySeeds.size))

		for iterSeed in xrange(mySeeds.size):
			seed = mySeeds.astype("int")[iterSeed]
			print "%d" % seed
			tc.randStream.seed = seed
			tc.rnaPartition.countsBulkIs(
				initRnaCnts * numpy.ones(tc.rnaPartition.countsBulk().shape)
				)

			ntpUsage = numpy.zeros((tc.metabolitePartition.ntps.countsBulk().size, lengthSec))
			rnaProduction = numpy.zeros((tc.rnaPartition.countsBulk().size, lengthSec))

			for t in xrange(lengthSec):
				ntps = ntpCounts * numpy.ones(tc.metabolitePartition.ntps.countsBulk().shape)
				tc.metabolitePartition.ntps.countsBulkIs(ntps)

				tc.enzymePartition.countsBulkIs(
					numpy.round(initEnzCnts * numpy.exp(numpy.log(2) / T_d * t))
					* numpy.ones(tc.enzymePartition.countsBulk().shape)
					)

				tc.evolveState()
				ntpUsage[:, t] = ntps - tc.metabolitePartition.ntps.countsBulk()
				rnaProduction[:, t] = tc.rnaPartition.countsBulk().flatten()

			myNtpUsage[:, :, iterSeed] = ntpUsage
			myRnaProduction[:, :, iterSeed] = rnaProduction

			fNtpUsage = ntpUsage / numpy.sum(ntpUsage, axis = 0)
			print "%d, %s" % (comm.rank, str(numpy.mean(fNtpUsage, axis = 1)))
			rnaMassProduction = numpy.diff(numpy.dot(rnaMWs, rnaProduction / Constants.nAvogadro))

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

			# Assert exponential increase in production of rRNA 23S
			rrlIdxs = tc.rnaPartition._getIndices(["RRLA-RRNA:nascent[c]",
				"RRLB-RRNA:nascent[c]", "RRLC-RRNA:nascent[c]",
				"RRLD-RRNA:nascent[c]", "RRLE-RRNA:nascent[c]",
				"RRLG-RRNA:nascent[c]", "RRLH-RRNA:nascent[c]"])[0]
			self.assertTrue(numpy.allclose(1., numpy.mean(nProd[rrlIdxs] / E[rrlIdxs]), rtol = 0, atol = 1e-1))

			# Loosely assert exponential increase in production all RNAs (this is noisy)
			numpyHandleError = numpy.geterr()
			numpy.seterr(divide = "ignore", invalid = "ignore")
			self.assertTrue(scipy.stats.nanmean(nProd / E) > 0.8)
			numpy.seterr(**numpyHandleError)

		## MPI communications ##
		# Get ready to gather values
		if comm.rank == 0:
			allNtpUsage = numpy.zeros(tc.metabolitePartition.ntps.countsBulk().size * lengthSec * allSeeds.size)
			allRnaProduction = numpy.zeros(tc.rnaPartition.countsBulk().size * lengthSec * allSeeds.size)

		else:
			allNtpUsage = None
			allRnaProduction = None

		# Gather NTP usage values
		recvcounts = sendcounts * ntpUsage.size
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
		comm.Gatherv(numpy.concatenate(numpy.split(myNtpUsage, mySeeds.size, axis = 2)), [allNtpUsage, recvcounts, displacements, MPI.DOUBLE])

		# Gather RNA production values
		recvcounts = sendcounts * rnaProduction.size
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
		comm.Gatherv(numpy.concatenate(numpy.split(myRnaProduction, mySeeds.size, axis = 2)), [allRnaProduction, recvcounts, displacements, MPI.DOUBLE])

		# Assemble gathered values
		if comm.rank == 0:
			allNtpUsage = numpy.dstack([x.reshape(ntpUsage.shape) for x in numpy.split(allNtpUsage.reshape(-1), allSeeds.size)])
			allRnaProduction = numpy.dstack([x.reshape(rnaProduction.shape) for x in numpy.split(allRnaProduction.reshape(-1), allSeeds.size)])
			
			# Check for proper assembly of the gathered data
			# (Check that it's not incorrect, not necessarily guaranteed to be correct)
			self.assertTrue(numpy.all(allNtpUsage[:, :, :mySeeds.size] == myNtpUsage))
			self.assertTrue(numpy.all(allRnaProduction[:, :, :mySeeds.size] == myRnaProduction))

	@noseAttrib.attr("rnaTotalProduction")
	def test_total_production(self):
		if comm.rank != 0:
			return

		sim = self.sim
		tc = sim.getProcess("Transcription")
		rnaMWs = tc.rnaPartition._state._massSingle.flat[tc.rnaPartition.mapping]

		rnaMWs[rnaMWs < 0] = 0

		ntpCounts = 1e6
		initEnzCnts = 960.
		initRnaCnts = 0.
		T_d = 3600.
		lengthSec = int(T_d)

		nSeeds = 100
		nProc = comm.size
		sendcounts = (nSeeds / nProc) * numpy.ones(nProc, dtype = int) + (numpy.arange(nProc, dtype = int) < nSeeds % nProc)
		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(sendcounts)[:-1]])
		mySeeds = numpy.zeros(sendcounts[comm.rank])
		allSeeds = numpy.arange(nSeeds, dtype = float)
		comm.Scatterv([allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE], mySeeds)
		print "Rank %d mySeeds: %s" % (comm.rank, str(mySeeds))
		
		# myNtpUsage = numpy.zeros((tc.metabolitePartition.idx["ntps"].size, lengthSec, mySeeds.size))
		# myRnaProduction = numpy.zeros((tc.rnaPartition.countsBulk().size, lengthSec, mySeeds.size))

		for iterSeed in xrange(mySeeds.size):
			seed = mySeeds.astype("int")[iterSeed]
			print "%d" % seed
			tc.randStream.seed = seed
			tc.rnaPartition.countsBulkIs(
				initRnaCnts * numpy.ones(tc.rnaPartition.countsBulk().shape)
				)

			ntpUsage = numpy.zeros((tc.metabolitePartition.ntps.countsBulk().size, lengthSec))
			rnaProduction = numpy.zeros((tc.rnaPartition.countsBulk().size, lengthSec))

			for t in xrange(lengthSec):
				ntps = ntpCounts * numpy.ones(tc.metabolitePartition.ntps.countsBulk().shape)
				tc.metabolitePartition.ntps.countsBulkIs(ntps)

				tc.enzymePartition.countsBulkIs(
					numpy.round(initEnzCnts * numpy.exp(numpy.log(2) / T_d * t))
					* numpy.ones(tc.enzymePartition.countsBulk().shape)
					)

				tc.evolveState()
				ntpUsage[:, t] = ntps - tc.metabolitePartition.ntps.countsBulk()
				rnaProduction[:, t] = tc.rnaPartition.countsBulk().flatten()

			# myNtpUsage[:, :, iterSeed] = ntpUsage
			# myRnaProduction[:, :, iterSeed] = rnaProduction

			fNtpUsage = ntpUsage / numpy.sum(ntpUsage, axis = 0)
			print "%d, %s" % (comm.rank, str(numpy.mean(fNtpUsage, axis = 1)))
			rnaMassProduction = numpy.diff(numpy.dot(rnaMWs, rnaProduction / Constants.nAvogadro))

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
			rrlIdxs = tc.rnaPartition._getIndices(["RRLA-RRNA:nascent[c]",
				"RRLB-RRNA:nascent[c]", "RRLC-RRNA:nascent[c]",
				"RRLD-RRNA:nascent[c]", "RRLE-RRNA:nascent[c]",
				"RRLG-RRNA:nascent[c]", "RRLH-RRNA:nascent[c]"])[0]
			self.assertTrue(numpy.allclose(1., numpy.mean(nProd[rrlIdxs] / E[rrlIdxs]), rtol = 0, atol = 1e-1))

			# Loosely assert exponential increase in all RNAs (this is noisy)
			numpyHandleError = numpy.geterr()
			numpy.seterr(divide = "ignore", invalid = "ignore")
			self.assertTrue(scipy.stats.nanmean(nProd / E) > 0.8)
			numpy.seterr(**numpyHandleError)
