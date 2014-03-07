"""
Test the interactions of Transcription.py, RnaDegradation.py, and Translation.py
Examines interplay between Transcription, Rna Degradation, and Translation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/26/2013
"""

from __future__ import division

# import unittest
# import warnings
# import nose.plugins.attrib as noseAttrib

# import numpy
# import scipy.stats
# import cPickle
# import os
# import copy
# #import matplotlib
# # matplotlib.use("agg")
# #from matplotlib import pyplot as plt
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

# class Test_Transcription_RnaDegradation(unittest.TestCase):

# 	@classmethod
# 	def setUpClass(cls):
# 		pass

# 	@classmethod
# 	def tearDownClass(cls):
# 		pass

# 	def setUp(self):
# 		self.sim = None
# 		self.outDir = "out/test/Test_Transcription_RnaDegradation_Translation"

# 		if comm.rank == 0:
# 			self.sim = cPickle.load(open(os.path.join("fixtures", "Simulation.cPickle"), "r"))
# 			if not os.path.exists(self.outDir):
# 				os.makedirs(self.outDir)

		
# 		self.sim = comm.bcast(self.sim, root = 0)	
# 		print "%s" % (self.sim.getState("Mass").meta["id"])

# 	def tearDown(self):
# 		pass


# 	# Tests
# 	@noseAttrib.attr('largetest')
# 	def test_net_rna_and_monomer_production(self):

# 		##### MPI Communications #####
# 		# Broadcast seeds
# 		nSeeds = 10
# 		nProc = comm.size
# 		sendcounts = (nSeeds / nProc) * numpy.ones(nProc, dtype = int) + (numpy.arange(nProc, dtype = int) < nSeeds % nProc)
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(sendcounts)[:-1]])
# 		mySeeds = numpy.zeros(sendcounts[comm.rank])
# 		allSeeds = numpy.arange(nSeeds, dtype = float)
# 		comm.Scatterv([allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE], mySeeds)

# 		sim = copy.deepcopy(self.sim)
# 		lengthSec = 3600

# 		tc = sim.getProcess("Transcription")
# 		rm = sim.getProcess("RnaMaturation")
# 		rd = sim.getProcess("RnaDegradation")
# 		tl = sim.getProcess("Translation")
# 		pm = sim.getProcess("ProteinMaturation")
# 		mc = sim.getState("MoleculeCounts")

# 		# Data to collect
# 		myTcRnaProd = numpy.zeros((tc.rna.counts.size, lengthSec, mySeeds.size))
# 		myRdRnaDegr = numpy.zeros((rd.rna.counts.size, lengthSec, mySeeds.size))
# 		myTcNtpUsage = numpy.zeros((tc.metabolite.parentState.tcNtpUsage.size, lengthSec, mySeeds.size))
# 		myRdNtpReturn = numpy.zeros((rd.metabolite.idx["ntps"].size, lengthSec, mySeeds.size))
# 		myTlMonProd = numpy.zeros((tl.protein.counts.size, lengthSec, mySeeds.size))
# 		myTlAaUsage = numpy.zeros((tl.metabolite.idx["aas"].size, lengthSec, mySeeds.size))
# 		myTotRnaCnts = numpy.zeros((mc.idx["matureRna"].size, lengthSec, mySeeds.size))
# 		myTotMonCnts = numpy.zeros((mc.idx["matureMonomers"].size, lengthSec, mySeeds.size))

# 		for iterSeed in xrange(mySeeds.size):

# 			seed = mySeeds.astype("int")[iterSeed]

# 			sim = copy.deepcopy(self.sim)
# 			sim.seed = seed
# 			sim.calcInitialConditions()

# 			tc = sim.getProcess("Transcription")
# 			rm = sim.getProcess("RnaMaturation")
# 			rd = sim.getProcess("RnaDegradation")
# 			tl = sim.getProcess("Translation")
# 			pm = sim.getProcess("ProteinMaturation")
# 			mc = sim.getState("MoleculeCounts")
# 			rd.rna.mws[rd.rna.mws < 0 ] = 0
# 			tc.rna.mws[tc.rna.mws < 0 ] = 0
# 			tl.protein.mws[tl.protein.mws < 0] = 0
# 			mc.mws[mc.mws < 0] = 0
			
# 			initRnapCnts = 1000.
# 			ntpCounts = 1e6

# 			initRnaseCnts = 1000.
# 			h2oCounts = 1e6

# 			initRibCnts = 8400
# 			aaCounts = 1e8

# 			T_d = 3600.
			

# 			initRnaMass = numpy.dot(mc.mws[mc.idx["matureRna"]], mc.counts[mc.idx["matureRna"], mc.cIdx["c"]]) / Constants.nAvogadro
# 			initRnaCnts = mc.counts[mc.idx["matureRna"], mc.cIdx["c"]].copy()
# 			initMonMass = numpy.dot(mc.mws[mc.idx["matureMonomers"]], numpy.sum(mc.counts[mc.idx["matureMonomers"], :], axis = 1)) / Constants.nAvogadro
# 			initMonCnts = numpy.sum(mc.counts[mc.idx["matureMonomers"], :], axis = 1)

# 			for t in xrange(lengthSec):
# 				mc.prepartition()
# 				mc.partition()

# 				tc.metabolite.counts[tc.metabolite.idx["ntps"]] = ntpCounts * numpy.ones(tc.metabolite.idx["ntps"].shape)
# 				tcEnz = tc.enzyme.counts.copy()
# 				tc.enzyme.counts = numpy.round(initRnapCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tc.enzyme.counts.shape)

# 				rd.metabolite.counts[rd.metabolite.idx["h2o"]] = h2oCounts
# 				rdEnz = rd.enzyme.counts.copy()
# 				rd.enzyme.counts = numpy.round(initRnaseCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(rd.enzyme.counts.shape)
# 				rdBegRnaCounts = rd.rna.counts.copy()

# 				tl.metabolite.counts[tl.metabolite.idx["aas"]] = aaCounts * numpy.ones(tl.metabolite.idx["aas"].shape)
# 				tlEnz = tl.enzyme.counts.copy()
# 				tl.enzyme.counts = numpy.round(initRibCnts / tl.enzyme.idx["23S"].size * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tl.enzyme.counts.shape)

# 				tc.evolveState()
# 				rm.evolveState()
# 				rd.evolveState()
# 				tl.evolveState()
# 				pm.evolveState()

# 				tc.enzyme.counts = tcEnz.copy()
# 				rd.enzyme.counts = rdEnz.copy()
# 				tl.enzyme.counts = tlEnz.copy()

# 				myTcRnaProd[:, t, iterSeed] = tc.rna.counts
# 				myRdRnaDegr[:, t, iterSeed] = rdBegRnaCounts - rd.rna.counts
# 				myTcNtpUsage[:, t, iterSeed] = tc.metabolite.parentState.tcNtpUsage
# 				myRdNtpReturn[:, t, iterSeed] = rd.metabolite.counts[rd.metabolite.idx["ntps"]]
# 				myTlMonProd[:, t, iterSeed] = tl.protein.counts
# 				myTlAaUsage[:, t, iterSeed] = tl.aasUsed

# 				mc.merge()
# 				mc.calculate()

# 				myTotRnaCnts[:, t, iterSeed] = mc.counts[mc.idx["matureRna"], mc.cIdx["c"]]
# 				myTotMonCnts[:, t, iterSeed] = numpy.sum(mc.counts[mc.idx["matureMonomers"], :], axis = 1)

# 				if comm.rank == 0:
# 					print "===== Time: %4d =====" % (t)
# 					print "Transcription RNAs produced: %d" % (int(numpy.sum(myTcRnaProd[:, t, iterSeed])))
# 					print "RnaDegradation RNAs degraded: %d" % (int(numpy.sum(myRdRnaDegr[:, t, iterSeed])))
# 					print "Translation monomers produced: %d" % (int(numpy.sum(myTlMonProd[:, t, iterSeed])))
# 					print "Total RNA mass (fg): %0.3f " % (1e15 * numpy.dot(mc.mws[mc.idx["matureRna"]], myTotRnaCnts[:, t, iterSeed]) / Constants.nAvogadro)
# 					print "Expected RNA mass (fg): %0.3f" % (1e15 * initRnaMass * numpy.exp(numpy.log(2) / T_d * t))
# 					print "Total monomer mass (fg): %0.3f" % (1e15 * numpy.dot(mc.mws[mc.idx["matureMonomers"]], myTotMonCnts[:, t, iterSeed]) / Constants.nAvogadro)
# 					print "Expected monomer mass (fg): %0.3f" % (1e15 * initMonMass * numpy.exp(numpy.log(2) / T_d * t))
# 					print "\n"


# 			# TODO: get rid of these two lines?
# 			finalRnaMass = numpy.dot(mc.mws[mc.idx["matureRna"]], myTotRnaCnts[:, -1, iterSeed]) / Constants.nAvogadro
# 			finalRnaCnts = mc.counts[mc.idx["matureRna"], mc.cIdx["c"]].copy()

# 		##### MPI Communications #####
# 		# Get ready to gather values
# 		if comm.rank == 0:
# 			allTcRnaProd = numpy.zeros(tc.rna.counts.size * lengthSec * allSeeds.size)
# 			allRdRnaDegr = numpy.zeros(rd.rna.counts.size * lengthSec * allSeeds.size)
# 			allTcNtpUsage = numpy.zeros(tc.metabolite.parentState.tcNtpUsage.size * lengthSec * allSeeds.size)
# 			allRdNtpReturn = numpy.zeros(rd.metabolite.idx["ntps"].size * lengthSec * allSeeds.size)
# 			allTlMonProd = numpy.zeros(tl.protein.counts.size * lengthSec * allSeeds.size)
# 			allTlAaUsage = numpy.zeros(tl.metabolite.idx["aas"].size * lengthSec * allSeeds.size)
# 			allTotRnaCnts = numpy.zeros(mc.idx["matureRna"].size * lengthSec * allSeeds.size)
# 			allTotMonCnts = numpy.zeros(mc.idx["matureMonomers"].size * lengthSec * allSeeds.size)
# 		else:
# 			allTcRnaProd = None
# 			allRdRnaDegr = None
# 			allTcNtpUsage = None
# 			allRdNtpReturn = None
# 			allTlMonProd = None
# 			allTlAaUsage = None
# 			allTotRnaCnts = None
# 			allTotMonCnts = None

# 		# Gather tcRnaProd values
# 		recvcounts = sendcounts * myTcRnaProd.shape[0] * myTcRnaProd.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myTcRnaProd, mySeeds.size, axis = 2)), [allTcRnaProd, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather rdRnaDegr values
# 		recvcounts = sendcounts * myRdRnaDegr.shape[0] * myRdRnaDegr.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myRdRnaDegr, mySeeds.size, axis = 2)), [allRdRnaDegr, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather tcNtpUsage values
# 		recvcounts = sendcounts * myTcNtpUsage.shape[0] * myTcNtpUsage.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myTcNtpUsage, mySeeds.size, axis = 2)), [allTcNtpUsage, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather rdNtpReturn values
# 		recvcounts = sendcounts * myRdNtpReturn.shape[0] * myRdNtpReturn.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myRdNtpReturn, mySeeds.size, axis = 2)), [allRdNtpReturn, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather tlMonProd values
# 		recvcounts = sendcounts * myTlMonProd.shape[0] * myTlMonProd.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myTlMonProd, mySeeds.size, axis = 2)), [allTlMonProd, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather tlAaUsage values
# 		recvcounts = sendcounts * myTlAaUsage.shape[0] * myTlAaUsage.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myTlAaUsage, mySeeds.size, axis = 2)), [allTlAaUsage, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather totRnaCnts values
# 		recvcounts = sendcounts * myTotRnaCnts.shape[0] * myTotRnaCnts.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myTotRnaCnts, mySeeds.size, axis = 2)), [allTotRnaCnts, recvcounts, displacements, MPI.DOUBLE])

# 		# Gather totMonCnts values
# 		recvcounts = sendcounts * myTotMonCnts.shape[0] * myTotMonCnts.shape[1]
# 		displacements = numpy.hstack([numpy.zeros(1), numpy.cumsum(recvcounts)[:-1]])
# 		comm.Gatherv(numpy.concatenate(numpy.split(myTotMonCnts, mySeeds.size, axis = 2)), [allTotMonCnts, recvcounts, displacements, MPI.DOUBLE])

# 		if comm.rank != 0:
# 			return

# 		allTcRnaProd = numpy.dstack([x.reshape(myTcRnaProd.shape) for x in numpy.split(allTcRnaProd.reshape(-1), allSeeds.size)])
# 		allRdRnaDegr = numpy.dstack([x.reshape(myRdRnaDegr.shape) for x in numpy.split(allRdRnaDegr.reshape(-1), allSeeds.size)])
# 		allTcNtpUsage = numpy.dstack([x.reshape(myTcNtpUsage.shape) for x in numpy.split(allTcNtpUsage.reshape(-1), allSeeds.size)])
# 		allRdNtpReturn = numpy.dstack([x.reshape(myRdNtpReturn.shape) for x in numpy.split(allRdNtpReturn.reshape(-1), allSeeds.size)])
# 		allTlMonProd = numpy.dstack([x.reshape(myTlMonProd.shape) for x in numpy.split(allTlMonProd.reshape(-1), allSeeds.size)])
# 		allTlAaUsage = numpy.dstack([x.reshape(myTlAaUsage.shape) for x in numpy.split(allTlAaUsage.reshape(-1), allSeeds.size)])
# 		allTotRnaCnts = numpy.dstack([x.reshape(myTotRnaCnts.shape) for x in numpy.split(allTotRnaCnts.reshape(-1), allSeeds.size)])
# 		allTotMonCnts = numpy.dstack([x.reshape(myTotMonCnts.shape) for x in numpy.split(allTotMonCnts.reshape(-1), allSeeds.size)])

# 		self.assertTrue(numpy.all(allTcRnaProd[:, :, :mySeeds.size] == myTcRnaProd))
# 		self.assertTrue(numpy.all(allRdRnaDegr[:, :, :mySeeds.size] == myRdRnaDegr))
# 		self.assertTrue(numpy.all(allTcNtpUsage[:, :, :mySeeds.size] == myTcNtpUsage))
# 		self.assertTrue(numpy.all(allRdNtpReturn[:, :, :mySeeds.size] == myRdNtpReturn))
# 		self.assertTrue(numpy.all(allTlMonProd[:, :, :mySeeds.size] == myTlMonProd))
# 		self.assertTrue(numpy.all(allTlAaUsage[:, :, :mySeeds.size] == myTlAaUsage))
# 		self.assertTrue(numpy.all(allTotRnaCnts[:, :, :mySeeds.size] == myTotRnaCnts))
# 		self.assertTrue(numpy.all(allTotMonCnts[:, :, :mySeeds.size] == myTotMonCnts))

# 		numpy.savez("/data/data.npz", allTcRnaProd = allTcRnaProd, allTlMonProd = allTlMonProd, allRdRnaDegr = allRdRnaDegr, allTcNtpUsage = allTcNtpUsage, allRdNtpReturn = allRdNtpReturn, allTlAaUsage = allTlAaUsage, allTotRnaCnts = allTotRnaCnts, allTotMonCnts = allTotMonCnts)

# 		# import ipdb
# 		# ipdb.set_trace()

# 		# Assert that, on average, all monomer counts have doubled
# 		self.assertTrue(numpy.allclose(2., numpy.mean(numpy.mean(allTotMonCnts[numpy.all(allTotMonCnts[:, 0, :] > 0, axis = 1), -1, :] / allTotMonCnts[numpy.all(allTotMonCnts[:, 0, :] > 0, axis = 1), 0, :], axis = 1)), atol = 8e-2, rtol = 0.))

# 		# Plot RNA mass fractions for each sim
# 		rnaIds = [x[1] for x in enumerate(mc.ids) if x[0] in mc.idx["matureRna"]]
# 		idx = {}
# 		idx["totRnaCnts"] = {}
# 		idx["totRnaCnts"]["rRna23Ss"] = mc.idx["rRna23Ss"] - mc.idx["matureRna"][0]
# 		idx["totRnaCnts"]["rRna16Ss"] = mc.idx["rRna16Ss"] - mc.idx["matureRna"][0]
# 		idx["totRnaCnts"]["rRna5Ss"] = mc.idx["rRna5Ss"] - mc.idx["matureRna"][0]
# 		idx["totRnaCnts"]["rRnas"] = numpy.concatenate((idx["totRnaCnts"]["rRna23Ss"], idx["totRnaCnts"]["rRna16Ss"], idx["totRnaCnts"]["rRna5Ss"]))
# 		idx["totRnaCnts"]["tRnas"] = mc.idx["tRnas"] - mc.idx["matureRna"][0]
# 		idx["totRnaCnts"]["mRnas"] = mc.idx["matureMrna"] - mc.idx["matureRna"][0]
# 		idx["mc"] = {}
# 		idx["mc"]["rRnas"] = numpy.concatenate((mc.idx["rRna23Ss"], mc.idx["rRna16Ss"], mc.idx["rRna5Ss"]))
# 		time = numpy.arange(lengthSec)

# 		for iDepth in xrange(allTotRnaCnts.shape[2]):
# 			thisTotRnaCnts = allTotRnaCnts[:, :, iDepth]

# 			plt.figure(1)
# 			plt.clf()

# 			plt.subplot(2, 2, 1)
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureRna"]], thisTotRnaCnts) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureRna"]], thisTotRnaCnts[:, 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("Total RNA Mass (fg)")
# 			plt.legend(loc = "lower right")

# 			plt.subplot(2, 2, 2)
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[idx["mc"]["rRnas"]], thisTotRnaCnts[idx["totRnaCnts"]["rRnas"], :]) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[idx["mc"]["rRnas"]], thisTotRnaCnts[idx["totRnaCnts"]["rRnas"], 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("rRNA Mass (fg)")
# 			plt.legend(loc = "lower right")

# 			plt.subplot(2, 2, 3)
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["tRnas"]], thisTotRnaCnts[idx["totRnaCnts"]["tRnas"], :]) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["tRnas"]], thisTotRnaCnts[idx["totRnaCnts"]["tRnas"], 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("tRNA Mass (fg)")
# 			plt.legend(loc = "lower right")

# 			plt.subplot(2, 2, 4)
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureMrna"]], thisTotRnaCnts[idx["totRnaCnts"]["mRnas"], :]) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureMrna"]], thisTotRnaCnts[idx["totRnaCnts"]["mRnas"], 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("mRNA Mass (fg)")
# 			plt.legend(loc = "lower right")

# 			for axes in plt.gcf().get_axes():
# 				yax = axes.yaxis
# 				yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 				xax = axes.xaxis
# 				xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 			plt.tight_layout()

# 			plt.savefig(os.path.join(self.outDir, "RNA Mass.%02d.pdf" % iDepth))

# 		# Plot NTP usage for each sim
# 		for iDepth in xrange(allTcNtpUsage.shape[2]):
# 			thisTcNtpUsage = allTcNtpUsage[:, :, iDepth]

# 			plt.figure(2)
# 			plt.clf()

# 			plt.subplot(4, 1, 1)
# 			plt.plot(time, thisTcNtpUsage[0, :], linewidth = 2, color = "k")
# 			plt.ylabel("ATPs")

# 			plt.subplot(4, 1, 2)
# 			plt.plot(time, thisTcNtpUsage[1, :], linewidth = 2, color = "k")
# 			plt.ylabel("CTPs")

# 			plt.subplot(4, 1, 3)
# 			plt.plot(time, thisTcNtpUsage[2, :], linewidth = 2, color = "k")
# 			plt.ylabel("GTPs")

# 			plt.subplot(4, 1, 4)
# 			plt.plot(time, thisTcNtpUsage[3, :], linewidth = 2, color = "k")
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("UTPs")

# 			for axes in plt.gcf().get_axes():
# 				yax = axes.yaxis
# 				yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 				xax = axes.xaxis
# 				xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 			plt.tight_layout()

# 			plt.savefig(os.path.join(self.outDir, "Transcript NTP Usages.%02d.pdf" % iDepth))

# 			# Plot relative NTP usage by Transcription
# 			plt.figure(3)
# 			plt.clf()

# 			plt.subplot(1, 1, 1)
# 			h = plt.plot(time, numpy.transpose(thisTcNtpUsage / numpy.sum(thisTcNtpUsage, axis = 0)), linewidth = 1)
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("Transcription Relative NTP Usages")
# 			plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

# 			for axes in plt.gcf().get_axes():
# 				yax = axes.yaxis
# 				yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 				xax = axes.xaxis
# 				xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 			plt.tight_layout()

# 			plt.savefig(os.path.join(self.outDir, "Transcript NTP Usages Relative.%02d.pdf" % iDepth))

# 			# Assert that transcription's NTP usage matches expectation
# 			normalize = lambda x: numpy.array(x).astype("float") / numpy.sum(x)
# 			self.assertTrue(numpy.allclose(1., numpy.mean((thisTcNtpUsage / numpy.sum(thisTcNtpUsage, axis = 0)), axis = 1) / normalize(numpy.sum(tc.rnaSynthProb.reshape(-1, 1) * tc.rnaNtCounts, axis = 0)), atol = 0, rtol = 1e-2))


# 		# Plot protein mass fractions for each sim
# 		for iDepth in xrange(allTotMonCnts.shape[2]):
# 			thisTotMonCnts = allTotMonCnts[:, :, iDepth]

# 			plt.figure(4)
# 			plt.clf()

# 			plt.subplot(1, 1, 1)
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureMonomers"]], thisTotMonCnts) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
# 			plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureMonomers"]], thisTotMonCnts[:, 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
# 			plt.xlabel("Time (s)")
# 			plt.ylabel("Total Monomer Mass (fg)")
# 			plt.legend(loc = "lower right")

# 			for axes in plt.gcf().get_axes():
# 				yax = axes.yaxis
# 				yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 				xax = axes.xaxis
# 				xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 			plt.savefig(os.path.join(self.outDir, "Monomer Mass.%02d.pdf" % iDepth))

# 		# Plot AA usage for each sim
# 		for iDepth in xrange(allTlAaUsage.shape[2]):
# 			thisTlAaUsage = allTlAaUsage[:, :, iDepth]

# 			aaDict = {	"ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4, "GLU": 5, "GLN": 6, "GLY": 7, "HIS": 8, "ILE": 9,
# 						"LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14, "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20}

# 			aaKeys = sorted(aaDict.keys())

# 			plt.figure(5)
# 			plt.clf()

# 			for iAA in xrange(len(aaKeys)):
# 				aa = aaKeys[iAA]
# 				plt.subplot(5, 4, iAA + 1)
# 				plt.plot(time, thisTlAaUsage[aaDict[aa], :], linewidth = 2, color = "k")
# 				plt.ylabel(aa)

# 			for axes in plt.gcf().get_axes():
# 				yax = axes.yaxis
# 				yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 				xax = axes.xaxis
# 				xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 			plt.savefig(os.path.join(self.outDir, "AA usage.%02d.pdf" % iDepth))

# 		# Plot some proteins which do double for each sim

# 		# Plot some proteins which don't double for each sim

		

# 		for iDepth in xrange(allTotMonCnts.shape[2]):
# 			thisTotMonCnts = numpy.squeeze(allTotMonCnts[:, :, iDepth])
			
# 		# Plot number of RNAs produced and degraded over time
# 		plt.figure(2)
# 		plt.clf()

# 		plt.subplot(2, 1, 1)
# 		plt.plot(time, numpy.sum(tcRnaProd, axis = 0), linewidth = 1, color = "k", label = "Actual")
# 		Nmvg = 100
# 		plt.plot(time[:-1 * Nmvg], numpy.convolve(numpy.sum(tcRnaProd, axis = 0), numpy.ones(Nmvg) / Nmvg)[(Nmvg - 1) : -1 * Nmvg], linewidth = 2, color = "0.75", label = "Moving Average [N = %d]" % Nmvg)
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("RNA production (counts)")
# 		plt.legend(loc = "upper left")


# 		plt.subplot(2, 1, 2)
# 		plt.plot(time, numpy.sum(rdRnaDegr, axis = 0), linewidth = 1, color = "k", label = "Actual")
# 		Nmvg = 100
# 		plt.plot(time[:-1 * Nmvg], numpy.convolve(numpy.sum(rdRnaDegr, axis = 0), numpy.ones(Nmvg) / Nmvg)[(Nmvg - 1) : -1 * Nmvg], linewidth = 2, color = "0.75", label = "Moving Average [N = %d]" % Nmvg)
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("RNA degradation (counts)")
# 		plt.legend(loc = "upper left")

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "RNA Production and Degradation.pdf"))

# 		# Plot production of 23S rRNA

# 		# Plot production of 16S rRNA

# 		# Plot production of 5S rRNA

# 		# Plot production of some tRNAs

# 		# Plot production of some mRNAs

# 		# Plot rRNAs
# 		plt.figure(3)
# 		plt.clf()

# 		plt.subplot(3, 1, 1)
# 		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["rRna23Ss"]]), linewidth = 1, color = "k")
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("23S rRNAs (counts)")

# 		plt.subplot(3, 1, 2)
# 		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["rRna16Ss"]]), linewidth = 1, color = "k")
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("16S rRNAs (counts)")

# 		plt.subplot(3, 1, 3)
# 		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["rRna5Ss"]]), linewidth = 1, color = "k")
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("5S rRNAs (counts)")

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "rRNA Counts.pdf"))

# 		# Plot total normalized counts of tRNAs
# 		plt.figure(4)
# 		plt.clf()

# 		plt.subplot(1, 1, 1)
# 		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["tRnas"]]) / totRnaCnts[idx["totRnaCnts"]["tRnas"], 0], linewidth = 1, color = "k")
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("tRNAs (normalized counts)")

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "tRNA Counts.pdf"))

# 		# Plot NTP usage by Transcription
# 		plt.figure(5)
# 		plt.clf()

# 		plt.subplot(4, 1, 1)
# 		plt.plot(time, tcNtpUsage[0, :], linewidth = 2, color = "k")
# 		plt.ylabel("ATPs")

# 		plt.subplot(4, 1, 2)
# 		plt.plot(time, tcNtpUsage[1, :], linewidth = 2, color = "k")
# 		plt.ylabel("CTPs")

# 		plt.subplot(4, 1, 3)
# 		plt.plot(time, tcNtpUsage[2, :], linewidth = 2, color = "k")
# 		plt.ylabel("GTPs")

# 		plt.subplot(4, 1, 4)
# 		plt.plot(time, tcNtpUsage[3, :], linewidth = 2, color = "k")
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("UTPs")

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "Transcript NTP Usages.pdf"))

# 		# Plot relative NTP usage by Transcription
# 		plt.figure(6)
# 		plt.clf()

# 		plt.subplot(1, 1, 1)
# 		h = plt.plot(time, numpy.transpose(tcNtpUsage / numpy.sum(tcNtpUsage, axis = 0)), linewidth = 1)
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("Transcription Relative NTP Usages")
# 		plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "Transcript NTP Usages Relative.pdf"))

# 		# Assert that transcription's NTP usage matches expectation
# 		normalize = lambda x: numpy.array(x).astype("float") / numpy.sum(x)
# 		self.assertTrue(numpy.allclose(1., numpy.mean((tcNtpUsage / numpy.sum(tcNtpUsage, axis = 0)), axis = 1) / normalize(numpy.sum(tc.rnaSynthProb.reshape(-1, 1) * tc.rnaNtCounts, axis = 0)), atol = 0, rtol = 1e-2))

# 		# Plot NTP (NMP someday) return by RnaDegradation
# 		plt.figure(7)
# 		plt.clf()

# 		plt.subplot(1, 1, 1)
# 		h = plt.plot(time, numpy.transpose(rdNtpReturn / numpy.sum(rdNtpReturn, axis = 0)), linewidth = 1)
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("RNA Degradation NTP Products Relative")
# 		plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "RNA Degradation NTP Products Relative.pdf"))


# 		# Plot NTP usage minus what's returned by RnaDegradation
# 		plt.figure(8)
# 		plt.clf()

# 		plt.subplot(1, 1, 1)
# 		metProd = tcNtpUsage - rdNtpReturn
# 		h = plt.plot(time, numpy.transpose(metProd / numpy.sum(metProd, axis = 0)), linewidth = 1)
# 		plt.xlabel("Time (s)")
# 		plt.ylabel("(Transcription NTP Usage - RNA Degradation Return) Relative")
# 		plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

# 		for axes in plt.gcf().get_axes():
# 			yax = axes.yaxis
# 			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
# 			xax = axes.xaxis
# 			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

# 		plt.tight_layout()

# 		plt.savefig(os.path.join(self.outDir, "Transcription NTP Usage - RNA Degradation Return Relative.pdf"))


# 		# Assert that translation's AA usage matches expectation
# 		# (except for selenocysteine)
# 		self.assertTrue(numpy.allclose(1., (numpy.mean((tlAaUsage / numpy.sum(tlAaUsage, axis = 0)), axis = 1) / normalize(numpy.dot(mc.monExp, tl.proteinAaCounts)))[[x for x in xrange(21) if x != 15]], atol = 0, rtol = 1e-2))

# 		# import ipdb
# 		# ipdb.set_trace()

# 		saveDict = {
# 			"initRnapCnts": initRnapCnts,
# 			"ntpCounts": ntpCounts,
# 			"initRnaseCnts": initRnaseCnts,
# 			"h2oCounts": h2oCounts,
# 			"T_d": T_d,
# 			"lengthSec": lengthSec,
# 			"initRnaMass": initRnaMass,
# 			"initRnaCnts": initRnaCnts,
# 			"tcRnaProd": tcRnaProd,
# 			"tlMonProd": tlMonProd,
# 			"rdRnaDegr": rdRnaDegr,
# 			"tcNtpUsage": tcNtpUsage,
# 			"rdNtpReturn": rdNtpReturn,
# 			"tlAaUsage": tlAaUsage,
# 			"totRnaCnts": totRnaCnts,
# 			"totMonCnts": totMonCnts,
# 		}

# 		cPickle.dump(saveDict, open(os.path.join(self.outDir, "data.cPickle"), "w"), cPickle.HIGHEST_PROTOCOL)
# 		cPickle.dump(saveDict, open(os.path.join("/data", "data.cPickle"), "w"), cPickle.HIGHEST_PROTOCOL)
# 		return
