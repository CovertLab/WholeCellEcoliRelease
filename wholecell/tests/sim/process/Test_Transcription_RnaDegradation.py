"""
Test the interactions of Transcription.py and RnaDegradation.py
Examines interplay between Transcription and Rna Degradation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/26/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

import numpy
import scipy.stats
import cPickle
import os
#import matplotlib
#matplotlib.use("agg")
#from matplotlib import pyplot as plt
from wholecell.util.Constants import Constants

from mpi4py import MPI

comm = MPI.COMM_WORLD

import copy_reg
import types

class Test_Transcription_RnaDegradation(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = None
		self.outDir = "out/test/Test_Transcription_RnaDegradation"

		if comm.rank == 0:
			self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
			if not os.path.exists(self.outDir):
				os.makedirs(self.outDir)

		
		self.sim = comm.bcast(self.sim, root = 0)	
		print "%s" % (self.sim.states["Mass"].meta["id"])

	def tearDown(self):
		pass


	# Tests
	@noseAttrib.attr('largetest')
	def test_net_rna_production(self):
		sim = self.sim
		tc = sim.processes["Transcription"]
		rm = sim.processes["RnaMaturation"]
		rd = sim.processes["RnaDegradation"]
		mc = sim.states["MoleculeCounts"]
		rd.rna.mws[rd.rna.mws < 0 ] = 0
		tc.rna.mws[tc.rna.mws < 0 ] = 0
		mc.mws[mc.mws < 0] = 0
		
		initRnapCnts = 1000.
		ntpCounts = 1e6

		initRnaseCnts = 1000.
		h2oCounts = 1e6

		T_d = 3600.
		lengthSec = 3600

		initRnaMass = numpy.dot(mc.mws[mc.idx["matureRna"]], mc.counts[mc.idx["matureRna"], mc.cIdx["c"]]) / Constants.nAvogadro
		initRnaCnts = mc.counts[mc.idx["matureRna"], mc.cIdx["c"]].copy()

		tcRnaProd = numpy.zeros((tc.rna.counts.size, lengthSec))
		rdRnaDegr = numpy.zeros((rd.rna.counts.size, lengthSec))
		tcNtpUsage = numpy.zeros((tc.metabolite.parentState.tcNtpUsage.size, lengthSec))
		rdNtpReturn = numpy.zeros((rd.metabolite.idx["ntps"].size, lengthSec))
		totRnaCnts = numpy.zeros((mc.idx["matureRna"].size, lengthSec))

		for t in xrange(lengthSec):
			mc.prepartition()
			mc.partition()

			tc.metabolite.counts[tc.metabolite.idx["ntps"]] = ntpCounts * numpy.ones(tc.metabolite.idx["ntps"].shape)
			tc.enzyme.counts = numpy.round(initRnapCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(tc.enzyme.counts.shape)

			rd.metabolite.counts[rd.metabolite.idx["h2o"]] = h2oCounts
			rd.enzyme.counts = numpy.round(initRnaseCnts * numpy.exp(numpy.log(2) / T_d * t)) * numpy.ones(rd.enzyme.counts.shape)
			rdBegRnaCounts = rd.rna.counts.copy()

			tc.evolveState()
			rm.evolveState()
			rd.evolveState()

			tcRnaProd[:, t] = tc.rna.counts
			rdRnaDegr[:, t] = rdBegRnaCounts - rd.rna.counts
			tcNtpUsage[:, t] = tc.metabolite.parentState.tcNtpUsage
			rdNtpReturn[:, t] = rd.metabolite.counts[rd.metabolite.idx["ntps"]]

			mc.merge()
			mc.calculate()

			totRnaCnts[:, t] = mc.counts[mc.idx["matureRna"], mc.cIdx["c"]]

			print "===== Time: %4d =====" % (t)
			print "Transcription RNAs produced: %d" % (int(numpy.sum(tcRnaProd[:, t])))
			print "RnaDegradation RNAs degraded: %d" % (int(numpy.sum(rdRnaDegr[:, t])))
			print "Total RNA mass (fg): %0.3f " % (1e15 * numpy.dot(mc.mws[mc.idx["matureRna"]], totRnaCnts[:, t]) / Constants.nAvogadro)
			print "Expected RNA mass (fg): %0.3f" % (1e15 * initRnaMass * numpy.exp(numpy.log(2) / T_d * t))
			print "\n"


		finalRnaMass = numpy.dot(mc.mws[mc.idx["matureRna"]], totRnaCnts[:, -1]) / Constants.nAvogadro
		finalRnaCnts = mc.counts[mc.idx["matureRna"], mc.cIdx["c"]].copy()

		rnaIds = [x[1] for x in enumerate(mc.ids) if x[0] in mc.idx["matureRna"]]
		idx = {}
		idx["totRnaCnts"] = {}
		idx["totRnaCnts"]["rRna23Ss"] = mc.idx["rRna23Ss"] - mc.idx["matureRna"][0]
		idx["totRnaCnts"]["rRna16Ss"] = mc.idx["rRna16Ss"] - mc.idx["matureRna"][0]
		idx["totRnaCnts"]["rRna5Ss"] = mc.idx["rRna5Ss"] - mc.idx["matureRna"][0]
		idx["totRnaCnts"]["rRnas"] = numpy.concatenate((idx["totRnaCnts"]["rRna23Ss"], idx["totRnaCnts"]["rRna16Ss"], idx["totRnaCnts"]["rRna5Ss"]))
		idx["totRnaCnts"]["tRnas"] = mc.idx["tRnas"] - mc.idx["matureRna"][0]
		idx["totRnaCnts"]["mRnas"] = mc.idx["matureMrna"] - mc.idx["matureRna"][0]
		idx["mc"] = {}
		idx["mc"]["rRnas"] = numpy.concatenate((mc.idx["rRna23Ss"], mc.idx["rRna16Ss"], mc.idx["rRna5Ss"]))
		time = numpy.arange(lengthSec)

		# Plot RNA masses over time (+exponential)
		plt.figure(1)
		plt.clf()

		plt.subplot(2, 2, 1)
		plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureRna"]], totRnaCnts) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
		plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureRna"]], totRnaCnts[:, 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
		plt.xlabel("Time (s)")
		plt.ylabel("Total RNA Mass (fg)")
		plt.legend(loc = "lower right")

		plt.subplot(2, 2, 2)
		plt.plot(time, 1e15 * numpy.dot(mc.mws[idx["mc"]["rRnas"]], totRnaCnts[idx["totRnaCnts"]["rRnas"], :]) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
		plt.plot(time, 1e15 * numpy.dot(mc.mws[idx["mc"]["rRnas"]], totRnaCnts[idx["totRnaCnts"]["rRnas"], 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
		plt.xlabel("Time (s)")
		plt.ylabel("rRNA Mass (fg)")
		plt.legend(loc = "lower right")

		plt.subplot(2, 2, 3)
		plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["tRnas"]], totRnaCnts[idx["totRnaCnts"]["tRnas"], :]) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
		plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["tRnas"]], totRnaCnts[idx["totRnaCnts"]["tRnas"], 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
		plt.xlabel("Time (s)")
		plt.ylabel("tRNA Mass (fg)")
		plt.legend(loc = "lower right")

		plt.subplot(2, 2, 4)
		plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureMrna"]], totRnaCnts[idx["totRnaCnts"]["mRnas"], :]) / Constants.nAvogadro, linewidth = 2, color = "k", label = "Actual")
		plt.plot(time, 1e15 * numpy.dot(mc.mws[mc.idx["matureMrna"]], totRnaCnts[idx["totRnaCnts"]["mRnas"], 0]) / Constants.nAvogadro * numpy.exp(numpy.log(2) / T_d * time), linestyle = "-.", color = "0.25", label = "Exponential")
		plt.xlabel("Time (s)")
		plt.ylabel("mRNA Mass (fg)")
		plt.legend(loc = "lower right")

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "RNA Mass.pdf"))

		# Plot number of RNAs produced and degraded over time
		plt.figure(2)
		plt.clf()

		plt.subplot(2, 1, 1)
		plt.plot(time, numpy.sum(tcRnaProd, axis = 0), linewidth = 1, color = "k", label = "Actual")
		Nmvg = 100
		plt.plot(time[:-1 * Nmvg], numpy.convolve(numpy.sum(tcRnaProd, axis = 0), numpy.ones(Nmvg) / Nmvg)[(Nmvg - 1) : -1 * Nmvg], linewidth = 2, color = "0.75", label = "Moving Average [N = %d]" % Nmvg)
		plt.xlabel("Time (s)")
		plt.ylabel("RNA production (counts)")
		plt.legend(loc = "upper left")


		plt.subplot(2, 1, 2)
		plt.plot(time, numpy.sum(rdRnaDegr, axis = 0), linewidth = 1, color = "k", label = "Actual")
		Nmvg = 100
		plt.plot(time[:-1 * Nmvg], numpy.convolve(numpy.sum(rdRnaDegr, axis = 0), numpy.ones(Nmvg) / Nmvg)[(Nmvg - 1) : -1 * Nmvg], linewidth = 2, color = "0.75", label = "Moving Average [N = %d]" % Nmvg)
		plt.xlabel("Time (s)")
		plt.ylabel("RNA degradation (counts)")
		plt.legend(loc = "upper left")

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "RNA Production and Degradation.pdf"))

		# Plot production of 23S rRNA

		# Plot production of 16S rRNA

		# Plot production of 5S rRNA

		# Plot production of some tRNAs

		# Plot production of some mRNAs

		# Plot rRNAs
		plt.figure(3)
		plt.clf()

		plt.subplot(3, 1, 1)
		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["rRna23Ss"]]), linewidth = 1, color = "k")
		plt.xlabel("Time (s)")
		plt.ylabel("23S rRNAs (counts)")

		plt.subplot(3, 1, 2)
		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["rRna16Ss"]]), linewidth = 1, color = "k")
		plt.xlabel("Time (s)")
		plt.ylabel("16S rRNAs (counts)")

		plt.subplot(3, 1, 3)
		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["rRna5Ss"]]), linewidth = 1, color = "k")
		plt.xlabel("Time (s)")
		plt.ylabel("5S rRNAs (counts)")

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "rRNA Counts.pdf"))

		# Plot total normalized counts of tRNAs
		plt.figure(4)
		plt.clf()

		plt.subplot(1, 1, 1)
		plt.plot(time, numpy.transpose(totRnaCnts[idx["totRnaCnts"]["tRnas"]]) / totRnaCnts[idx["totRnaCnts"]["tRnas"], 0], linewidth = 1, color = "k")
		plt.xlabel("Time (s)")
		plt.ylabel("tRNAs (normalized counts)")

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "tRNA Counts.pdf"))

		# Plot NTP usage by Transcription
		plt.figure(5)
		plt.clf()

		plt.subplot(4, 1, 1)
		plt.plot(time, tcNtpUsage[0, :], linewidth = 2, color = "k")
		plt.ylabel("ATPs")

		plt.subplot(4, 1, 2)
		plt.plot(time, tcNtpUsage[1, :], linewidth = 2, color = "k")
		plt.ylabel("CTPs")

		plt.subplot(4, 1, 3)
		plt.plot(time, tcNtpUsage[2, :], linewidth = 2, color = "k")
		plt.ylabel("GTPs")

		plt.subplot(4, 1, 4)
		plt.plot(time, tcNtpUsage[3, :], linewidth = 2, color = "k")
		plt.xlabel("Time (s)")
		plt.ylabel("UTPs")

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "Transcript NTP Usages.pdf"))

		# Plot relative NTP usage by Transcription
		plt.figure(6)
		plt.clf()

		plt.subplot(1, 1, 1)
		h = plt.plot(time, numpy.transpose(tcNtpUsage / numpy.sum(tcNtpUsage, axis = 0)), linewidth = 1)
		plt.xlabel("Time (s)")
		plt.ylabel("Transcription Relative NTP Usages")
		plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "Transcript NTP Usages Relative.pdf"))

		# Assert that transcription's NTP usage matches expectation
		normalize = lambda x: numpy.array(x).astype("float") / numpy.sum(x)
		self.assertTrue(numpy.allclose(1., numpy.mean((tcNtpUsage / numpy.sum(tcNtpUsage, axis = 0)), axis = 1) / normalize(numpy.sum(tc.rnaSynthProb.reshape(-1, 1) * tc.rnaNtCounts, axis = 0)), atol = 0, rtol = 1e-2))

		# Plot NTP (NMP someday) return by RnaDegradation
		plt.figure(7)
		plt.clf()

		plt.subplot(1, 1, 1)
		h = plt.plot(time, numpy.transpose(rdNtpReturn / numpy.sum(rdNtpReturn, axis = 0)), linewidth = 1)
		plt.xlabel("Time (s)")
		plt.ylabel("RNA Degradation NTP Products Relative")
		plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "RNA Degradation NTP Products Relative.pdf"))


		# Plot NTP usage minus what's returned by RnaDegradation
		plt.figure(8)
		plt.clf()

		plt.subplot(1, 1, 1)
		metProd = tcNtpUsage - rdNtpReturn
		h = plt.plot(time, numpy.transpose(metProd / numpy.sum(metProd, axis = 0)), linewidth = 1)
		plt.xlabel("Time (s)")
		plt.ylabel("(Transcription NTP Usage - RNA Degradation Return) Relative")
		plt.legend(h, ["ATP", "CTP", "GTP", "UTP"])

		for axes in plt.gcf().get_axes():
			yax = axes.yaxis
			yax.set_ticks([yax.get_ticklocs()[0], yax.get_ticklocs()[-1]])
			xax = axes.xaxis
			xax.set_ticks([xax.get_ticklocs()[0], xax.get_ticklocs()[-1]])

		plt.tight_layout()

		plt.savefig(os.path.join(self.outDir, "Transcription NTP Usage - RNA Degradation Return Relative.pdf"))

		# import ipdb
		# ipdb.set_trace()

		saveDict = {
			"initRnapCnts": initRnapCnts,
			"ntpCounts": ntpCounts,
			"initRnaseCnts": initRnaseCnts,
			"h2oCounts": h2oCounts,
			"T_d": T_d,
			"lengthSec": lengthSec,
			"initRnaMass": initRnaMass,
			"initRnaCnts": initRnaCnts,
			"tcRnaProd": tcRnaProd,
			"rdRnaDegr": rdRnaDegr,
			"tcNtpUsage": tcNtpUsage,
			"rdNtpReturn": rdNtpReturn,
		}

		cPickle.dump(saveDict, open(os.path.join(self.outDir, "data.cPickle"), "w"), cPickle.HIGHEST_PROTOCOL)

		return