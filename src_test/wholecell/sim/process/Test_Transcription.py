"""
Test Transcription.py
Tests Transcription process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/16/2013
"""

import unittest
import os

import nose.plugins.attrib as noseAttrib
import numpy
import tables
import scipy.stats

from wholecell.util.Constants import Constants

FIXTURE_DIR = os.path.join('out', 'test', 'rnaProduction')

# NOTE: see generateFitterTestFixtures.py for the simulations used in this file

class Test_Transcription(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		doublingTime = 3600.
		simLength = 500

		dirs = os.walk(FIXTURE_DIR).next()[1] # get all of the directories

		nSims = len(dirs)

		ntpIDs = ['ATP', 'UTP', 'CTP', 'GTP']
		rrnaIDs = ['RRLA-RRNA:nascent', 'RRLB-RRNA:nascent',
			'RRLC-RRNA:nascent', 'RRLD-RRNA:nascent', 'RRLE-RRNA:nascent',
			'RRLG-RRNA:nascent', 'RRLH-RRNA:nascent']

		assignedIdxs = False

		ntpIdxs = None
		rnaIdxs = None
		rrnaIdxs = None
		processIdx = None
		compartmentIdx = None

		width = 10

		cls.ntpInitialUsage = numpy.zeros((len(ntpIDs), nSims), float)
		cls.ntpFinalUsage = numpy.zeros((len(ntpIDs), nSims), float)

		cls.rnaInitialProduction = None
		cls.rnaFinalProduction = None

		cls.expectedRatio = numpy.exp(numpy.log(2)/doublingTime * simLength)

		for iSim, simDir in enumerate(dirs):
			with tables.openFile(os.path.join(FIXTURE_DIR, simDir, 'MoleculeCounts.hdf')) as h5file:
				if not assignedIdxs:
					names = h5file.get_node('/names')
					indexes = h5file.get_node('/indexes')

					molIDs = names.molIDs.read()
					processes = names.processes.read()
					compartments = names.compartments.read()
					
					ntpIdxs = numpy.array([molIDs.index(id_) for id_ in ntpIDs])
					rnaIdxs = indexes.nascentRnas.read()
					rrnaIdxs = numpy.array([molIDs.index(id_) for id_ in rrnaIDs])

					processIdx = processes.index('Transcription')
					compartmentIdx = compartments.index('c')

					cls.rnaInitialProduction = numpy.zeros((len(rnaIdxs), nSims), float)
					cls.rnaFinalProduction = numpy.zeros((len(rnaIdxs), nSims), float)

					cls.rnaFinalCount = numpy.zeros((len(rnaIdxs), nSims), float)

					assignedIdxs = True

				mc = h5file.root.MoleculeCounts

				cls.ntpInitialUsage[:, iSim] = (mc.read(1, 1+width, None, 'countsBulkPartitioned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(1, width+1, None, 'countsBulkReturned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0))
				cls.ntpFinalUsage[:, iSim] = (mc.read(simLength+1 - width, simLength+1, None, 'countsBulkPartitioned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(simLength+1 - width, simLength+1, None, 'countsBulkReturned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0))

				cls.rnaInitialProduction[:, iSim] = (mc.read(1, 1+width, None, 'countsBulkReturned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(1, width+1, None, 'countsBulkPartitioned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0))
				cls.rnaFinalProduction[:, iSim] = (mc.read(simLength+1 - width, simLength+1, None, 'countsBulkReturned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(simLength+1 - width, simLength+1, None, 'countsBulkPartitioned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0))

				cls.rnaFinalCount[:, iSim] = mc[-1]['countsBulk'][rnaIdxs, compartmentIdx]

				# indexing order for mc.read(...): time, molecule, compartment, partition

		assignedFittedParameters = False

		cls.rnaSynthProb = None

		for iSim, simDir in enumerate(dirs):
			with tables.openFile(os.path.join(FIXTURE_DIR, simDir, 'Main.hdf')) as h5file:
				if not assignedFittedParameters:
					cls.rnaSynthProb = h5file.get_node('/fitParameters').rnaSynthProb.read()

					assignedFittedParameters = True


	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		pass


	def tearDown(self):
		pass


	@noseAttrib.attr('largetest', 'modelfitting', 'transcription')
	def test_ntpUsage(self):
		# Test for exponential NTP usage
		# NOTE: a proper fit should actually exceed the expected ratio since this doesn't include degradation
		ntpRatio = self.ntpFinalUsage/self.ntpInitialUsage
		self.assertTrue(numpy.allclose(self.expectedRatio, ntpRatio, rtol = 0.25))
		self.assertTrue((self.expectedRatio <= ntpRatio).all(), 'Final NTP usage was less than expected.')


	@noseAttrib.attr('largetest', 'modelfitting', 'transcription')
	def test_rnaProduction(self):
		# Test for exponential RNA production
		# NOTE: a proper fit should actually exceed the expected ratio since this doesn't include degradation
		rnaRatio = self.rnaFinalProduction.sum()/self.rnaInitialProduction.sum()
		self.assertTrue(numpy.allclose(self.expectedRatio, rnaRatio, rtol = 0.25))
		self.assertTrue((self.expectedRatio <= rnaRatio).all(), 'Final RNA production was less than expected.')

		# Test for appropriate ratios of rRNA production
		averageCounts = self.rnaFinalCount.mean(1)

		expectedProduction = self.rnaSynthProb * averageCounts.sum()

		# Assert at least a 95% positive correlation
		self.assertGreater(
			scipy.stats.pearsonr(averageCounts, expectedProduction)[0],
			0.95
			)


	# Tests
	# @noseAttrib.attr('rnaProduction')
	# @noseAttrib.attr('largetest')
	def test_production_old(self):
		sim = self.sim
		tc = sim.processes["Transcription"]
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

	# @noseAttrib.attr("rnaTotalProduction")
	# @noseAttrib.attr('largetest')
	def test_total_production(self):
		if comm.rank != 0:
			return

		sim = self.sim
		tc = sim.processes["Transcription"]
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
