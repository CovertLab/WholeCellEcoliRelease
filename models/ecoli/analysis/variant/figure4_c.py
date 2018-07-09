"""
Plot to compare cell properties across different growth conditions similar to Dennis and Bremer. 1996. Fig 2
Multiple generations of the same variant will be plotted together

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/7/16
"""

from __future__ import absolute_import


import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

NT_MW = 487.0
PROTEIN_MW = 110.0

FONT_SIZE=9

def mm2inch(value):
	return value * 0.0393701

trim = 0.1

#size1 = (14.429, 18.9)
size1 = (29.25, 24.401)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, "inputDir does not currently exist as a directory"

		ap = AnalysisPaths(inputDir, variant_plot = True)
		all_cells = ap.get_cells(generation=[0,1], seed=[0])

		if ap.n_variant == 1:
			print "Disabled. Needs correct variant."
			return

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		rnaToProteinDict = {}
		dnaToProteinDict = {}
		elngRateDict = {}
		stableRnaFractionDict = {}
		mRnaFractionDict = {}
		doublingPerHourDict = {}
		stableRnaSynthRateDict = {}
		mRnaSynthRateDict = {}
		numOriginsAtInitDict = {}

		variantSimDataFile = ap.get_variant_kb(all_cells[0])
		sim_data = cPickle.load(open(variantSimDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro.asNumber()
		chromMass = (sim_data.getter.getMass(['CHROM_FULL[c]'])[0] / sim_data.constants.nAvogadro).asNumber()
		for simDir in all_cells:
			simOutDir = os.path.join(simDir, "simOut")
			variant = int(simDir[simDir.rfind('generation_')-14:simDir.rfind('generation_')-8])

			#print variant
			print "Loading sim"

			mass = TableReader(os.path.join(simOutDir, "Mass"))

			protein = mass.readColumn("proteinMass") * 10**-15
			rna = mass.readColumn("rnaMass") * 10**-15
			dna = mass.readColumn("dnaMass") * 10**-15

			growthRate = mass.readColumn("instantaniousGrowthRate")
			doublingTime = np.nanmean(np.log(2) / growthRate / 60)

			rnaNT = rna / NT_MW * nAvogadro
			proteinAA = protein / PROTEIN_MW * nAvogadro

			# Count chromosome equivalents
			chromEquivalents = dna / chromMass

			# Load ribosome data
			ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
			actualElongations = ribosomeDataFile.readColumn("actualElongations")
			ribosomeDataFile.close()

			transcriptDataFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
			rnaSynth = transcriptDataFile.readColumn("countRnaSynthesized")
			isTRna = sim_data.process.transcription.rnaData["isTRna"]
			isRRna = sim_data.process.transcription.rnaData["isRRna"]
			isMRna = sim_data.process.transcription.rnaData["isMRna"]
			stableRnaSynth = np.sum(rnaSynth[:,isTRna], axis=1) + np.sum(rnaSynth[:,isRRna], axis=1)
			totalRnaSynth = np.sum(rnaSynth, axis=1).astype(float)
			rnaFraction = stableRnaSynth / totalRnaSynth
			mRnaSynth = np.sum(rnaSynth[:,isMRna], axis=1)
			mRnaFraction = mRnaSynth/ totalRnaSynth

			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]


			originIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("originOfReplication")
			originOfReplication = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, originIndex]

			uniqueMoleculeCounts.close()

			criticalInitMass = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")

			numOriginsAtInit = originOfReplication[np.where(criticalInitMass >= 1.0)[0] - 1][0]


			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
			timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

			#import ipdb
			#ipdb.set_trace()

			if variant in rnaToProteinDict.keys():
				rnaToProteinDict[variant] = np.append(rnaToProteinDict[variant], rnaNT / (proteinAA / 100))
				dnaToProteinDict[variant] = np.append(dnaToProteinDict[variant], chromEquivalents / (proteinAA / 10**9))
				elngRateDict[variant] = np.append(elngRateDict[variant], (actualElongations / activeRibosome / timeStepSec)[3:])
				stableRnaFractionDict[variant] = np.append(stableRnaFractionDict[variant], np.asarray(rnaFraction)[~np.isnan(rnaFraction)])
				mRnaFractionDict[variant] = np.append(mRnaFractionDict[variant], np.asarray(mRnaFraction)[~np.isnan(mRnaFraction)])
				doublingPerHourDict[variant] = np.append(doublingPerHourDict[variant], 60 / doublingTime)
				stableRnaSynthRateDict[variant] = np.append(stableRnaSynthRateDict[variant], (np.asarray(stableRnaSynth)/timeStepSec/60/proteinAA/10**3))
				numOriginsAtInitDict[variant] = np.append(numOriginsAtInitDict[variant], numOriginsAtInit)
				#mRnaSynthRateDict[variant] = np.append(mRnaSynthRateDict[variant]), np.asarray(mRnaFraction)[~np.isnan(mRnaFraction)]/
			else:
				rnaToProteinDict[variant] = rnaNT / (proteinAA / 100)
				dnaToProteinDict[variant] = chromEquivalents / (proteinAA / 10**9)
				elngRateDict[variant] = (actualElongations / activeRibosome / timeStepSec)[3:]
				stableRnaFractionDict[variant] = np.asarray(rnaFraction)[~np.isnan(rnaFraction)]
				mRnaFractionDict[variant] = np.asarray(mRnaFraction)[~np.isnan(mRnaFraction)]
				doublingPerHourDict[variant] = 60 / doublingTime
				stableRnaSynthRateDict[variant] = (np.asarray(stableRnaSynth)/proteinAA/10**3)
				numOriginsAtInitDict[variant] = numOriginsAtInit

		#import ipdb
		#ipdb.set_trace()


		rnaToProtein = np.zeros(3)
		dnaToProtein = np.zeros(3)
		elngRate = np.zeros(3)
		stableRnaFraction = np.zeros(3)
		doublingPerHour = np.zeros(3)
		mRnaFractionCalc = np.zeros(3)
		rSSynthesisRate = np.zeros(3)
		stableRnaSynthRate = np.zeros(3)
		numOriginsAtInit = np.zeros(3)

		rnaToProtein_error = np.zeros(3)
		dnaToProtein_error = np.zeros(3)
		elngRate_error = np.zeros(3)
		stableRnaFraction_error = np.zeros(3)
		doublingPerHour_error = np.zeros(3)
		mRnaFractionCalc_error = np.zeros(3)
		rSSynthesisRate_error = np.zeros(3)
		stableRnaSynthRate_error = np.zeros(3)
		numOriginsAtInit_error = np.zeros(3)

		#get Averages:
		newOrder = [1,0,2]
		for j in range(0,3):
			rnaToProtein[j] = np.mean(rnaToProteinDict[newOrder[j]])
			dnaToProtein[j] = np.mean(dnaToProteinDict[newOrder[j]])
			elngRate[j] = np.mean(elngRateDict[newOrder[j]])
			stableRnaFraction[j] = np.mean(stableRnaFractionDict[newOrder[j]])
			mRnaFractionCalc[j] = np.mean(mRnaFractionDict[newOrder[j]])
			stableRnaSynthRate[j] = np.mean(stableRnaSynthRateDict[newOrder[j]])
			doublingPerHour[j] = np.mean(doublingPerHourDict[newOrder[j]])
			numOriginsAtInit[j] = np.mean(numOriginsAtInitDict[newOrder[j]])

			rnaToProtein_error[j] = np.std(rnaToProteinDict[newOrder[j]])
			dnaToProtein_error[j] = np.std(dnaToProteinDict[newOrder[j]])
			elngRate_error[j] = np.std(elngRateDict[newOrder[j]])
			stableRnaFraction_error[j] = np.std(stableRnaFractionDict[newOrder[j]])
			mRnaFractionCalc_error[j] = np.std(mRnaFractionDict[newOrder[j]])
			stableRnaSynthRate_error[j] = np.std(stableRnaSynthRateDict[newOrder[j]])
			doublingPerHour_error[j] = np.std(doublingPerHourDict[newOrder[j]])
			numOriginsAtInit_error[j] = np.std(numOriginsAtInitDict[newOrder[j]])

		#import all Dennis and Bremer Data:
		#dt order 0.6, 1.5, 2.5
		rnaToProteinDB = np.array([6.615, 10.962, 15.6])
		rnaToProteinDB_error = rnaToProteinDB * (0.06 + 0.025)
		dnaToProteinDB = np.array([2.82, 1.79, 1.52])
		dnaToProteinDB_error = dnaToProteinDB * (0.06 + 0.05)
		elngRateDB = np.array([12, 18, 21])
		elngRateDB_error = elngRateDB * (0.06 + 0.025)
		stableRnaFractionDB = np.array([.41,.68,.85])
		stableRnaFractionDB_error = stableRnaFractionDB * (0.06 + 0.025)
		mRnaFractionDB = np.array([.59, .32, .15])
		mRnaFractionDB_error = mRnaFractionDB * (0.06 + 0.025)
		stableRnaSynthRateDB = [0.53, 2.23, 5.3]
		numOriginsAtInitDB = np.array([1,2,4])
		numOriginsAtInitDB_error = numOriginsAtInitDB * (0.05 + 0.05)
		doublingPerHourDB = [0.6, 1.5, 2.5]


		############## PLOT 0 ###############

		mult = 3
		fig, ax0 = plt.subplots(1,1, figsize = (mm2inch(size1[0])*mult, mm2inch(size1[1])*mult))

		ax0.plot(doublingPerHour, rnaToProtein, color = 'k', linewidth = 1)
		ax0.errorbar(doublingPerHour, rnaToProtein, yerr=rnaToProtein_error, color = "k", linewidth = 1)
		ax0.plot(doublingPerHourDB, rnaToProteinDB, color = 'blue', linewidth = 1)
		ax0.errorbar(doublingPerHourDB, rnaToProteinDB, yerr = rnaToProteinDB_error, color = "blue", linewidth = 1)
		ax0.set_ylabel("RNA / Protein (nuc/100 aa)", fontsize=FONT_SIZE)
		ax0.set_ylim([0,20])

		whitePadSparklineAxis(ax0, False)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

		exportFigure(plt, plotOutDir, plotOutFileName + "0", metadata)

		for axes in [ax0]:
			axes.tick_params(
				axis='x',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				bottom='off',      # ticks along the bottom edge are off
				top='off',         # ticks along the top edge are off
				labelbottom='off') # labels along the bottom edge are off
			axes.tick_params(
				axis='y',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				left='off',      # ticks along the bottom edge are off
				right='off',         # ticks along the top edge are off
				labelleft='off') # labels along the bottom edge are off

			axes.set_xlabel("")
			axes.set_ylabel("")

		plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

		exportFigure(plt, plotOutDir, plotOutFileName + "0_stripped" ,metadata, transparent = True)
		plt.close("all")

		############## PLOT 1 ###############

		mult = 3
		fig, ax2 = plt.subplots(1,1, figsize = (mm2inch(size1[0])*mult, mm2inch(size1[1])*mult))

		ax2.plot(doublingPerHour, elngRate, color = 'k', linewidth = 1)
		ax2.errorbar(doublingPerHour, elngRate, yerr=elngRate_error, color = "k", linewidth = 1)
		ax2.plot(doublingPerHourDB, elngRateDB, color = 'blue', linewidth = 1)
		ax2.errorbar(doublingPerHourDB, elngRateDB, yerr = elngRateDB_error, color = "blue", linewidth = 1)
		ax2.set_ylabel("Ribosome Elongation Rate (aa/s)", fontsize=FONT_SIZE)
		ax2.set_ylim([0,25])

		whitePadSparklineAxis(ax2, False)

		for tick in ax2.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax2.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

		exportFigure(plt, plotOutDir, plotOutFileName + "1", metadata)

		for axes in [ax2]:
			axes.tick_params(
				axis='x',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				bottom='off',      # ticks along the bottom edge are off
				top='off',         # ticks along the top edge are off
				labelbottom='off') # labels along the bottom edge are off
			axes.tick_params(
				axis='y',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				left='off',      # ticks along the bottom edge are off
				right='off',         # ticks along the top edge are off
				labelleft='off') # labels along the bottom edge are off

			axes.set_xlabel("")
			axes.set_ylabel("")

		plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

		exportFigure(plt, plotOutDir, plotOutFileName + "1_stripped" ,metadata, transparent = True)
		plt.close("all")


		# ############## PLOT 2 ###############

		mult = 3
		fig, ax5 = plt.subplots(1,1, figsize = (mm2inch(size1[0])*mult, mm2inch(size1[1])*mult))

		ax5.plot(doublingPerHour,numOriginsAtInit, color='k', linewidth=1, label = 'Simulation')
		ax5.errorbar(doublingPerHour, numOriginsAtInit, yerr=numOriginsAtInit_error, color = "k", linewidth = 1)
		ax5.plot(doublingPerHourDB,numOriginsAtInitDB, color='blue', linewidth=1, label = 'Dennis and Bremer')
		ax5.errorbar(doublingPerHourDB,numOriginsAtInitDB, yerr=numOriginsAtInitDB_error, linewidth=1, color="blue")
		ax5.set_ylabel("Origins / Cell at Initiation", fontsize=FONT_SIZE)
		ax5.set_ylim([0,5])
		ax5.legend(loc=4,frameon=False, fontsize=FONT_SIZE)
		ax5.set_xlabel("Doublings per Hour", fontsize=FONT_SIZE)

		whitePadSparklineAxis(ax5, False)

		for tick in ax5.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax5.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

		exportFigure(plt, plotOutDir, plotOutFileName + "2", metadata)

		for axes in [ax5]:
			axes.tick_params(
				axis='x',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				bottom='off',      # ticks along the bottom edge are off
				top='off',         # ticks along the top edge are off
				labelbottom='off') # labels along the bottom edge are off
			axes.tick_params(
				axis='y',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				left='off',      # ticks along the bottom edge are off
				right='off',         # ticks along the top edge are off
				labelleft='off') # labels along the bottom edge are off

			axes.set_xlabel("")
			axes.set_ylabel("")

		plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

		exportFigure(plt, plotOutDir, plotOutFileName + "2_stripped" ,metadata, transparent = True)
		plt.close("all")

		# ############## PLOT 3 ###############

		mult = 3
		fig, ax3 = plt.subplots(1,1, figsize = (mm2inch(size1[0])*mult, mm2inch(size1[1])*mult))

		ax3.plot(doublingPerHour, stableRnaFraction, color = 'k', linewidth = 1)
		ax3.errorbar(doublingPerHour, stableRnaFraction, yerr=stableRnaFraction_error, color = "k", linewidth = 1)
		ax3.plot(doublingPerHourDB, stableRnaFractionDB, color = 'blue', linewidth = 1)
		ax3.errorbar(doublingPerHourDB, stableRnaFractionDB, yerr = stableRnaFractionDB_error, color = "blue", linewidth=1)
		ax3.set_ylabel("Synthesis rate Stable RNA / Total RNA")
		ax3.set_ylim([0,1])

		whitePadSparklineAxis(ax3, False)

		for tick in ax3.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax3.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

		exportFigure(plt, plotOutDir, plotOutFileName + "3", metadata)

		for axes in [ax3]:
			axes.tick_params(
				axis='x',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				bottom='off',      # ticks along the bottom edge are off
				top='off',         # ticks along the top edge are off
				labelbottom='off') # labels along the bottom edge are off
			axes.tick_params(
				axis='y',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				left='off',      # ticks along the bottom edge are off
				right='off',         # ticks along the top edge are off
				labelleft='off') # labels along the bottom edge are off

			axes.set_xlabel("")
			axes.set_ylabel("")

		plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

		exportFigure(plt, plotOutDir, plotOutFileName + "3_stripped" ,metadata, transparent = True)
		plt.close("all")





		# mult = 3
		# fig, ax = plt.subplots(3,1, sharex=True, figsize = (mm2inch(18)*mult, mm2inch(58)*mult))

		# ax0 = ax[0]
		# #ax1 = ax[0,1]
		# ax2 = ax[1]
		# #ax3 = ax[1,1]
		# #ax4 = ax[2,0]
		# ax5 = ax[2]

		# ax0.plot(doublingPerHour, rnaToProtein, color = 'k', linewidth = 1)
		# ax0.errorbar(doublingPerHour, rnaToProtein, yerr=rnaToProtein_error, color = "k", linewidth = 1)
		# ax0.plot(doublingPerHourDB, rnaToProteinDB, color = 'blue', linewidth = 1)
		# ax0.errorbar(doublingPerHourDB, rnaToProteinDB, yerr = rnaToProteinDB_error, color = "blue", linewidth = 1)
		# ax0.set_ylabel("RNA / Protein (nuc/100 aa)", fontsize=FONT_SIZE)
		# ax0.set_ylim([0,20])

		# # ax1.plot(doublingPerHour, dnaToProtein, color = 'k', linewidth = 1)
		# # ax1.errorbar(doublingPerHour, dnaToProtein, yerr=dnaToProtein_error, color = "k", linewidth = 1)
		# # ax1.plot(doublingPerHourDB, dnaToProteinDB, color = 'blue', linewidth = 1)
		# # ax1.errorbar(doublingPerHourDB, dnaToProteinDB, yerr=dnaToProteinDB_error, color="blue", linewidth = 1)
		# # ax1.set_ylabel("DNA / Protein (chrom eq/" + r"$10^9$" + " aa)")
		# # ax1.set_ylim([0,4])

		# ax2.plot(doublingPerHour, elngRate, color = 'k', linewidth = 1)
		# ax2.errorbar(doublingPerHour, elngRate, yerr=elngRate_error, color = "k", linewidth = 1)
		# ax2.plot(doublingPerHourDB, elngRateDB, color = 'blue', linewidth = 1)
		# ax2.errorbar(doublingPerHourDB, elngRateDB, yerr = elngRateDB_error, color = "blue", linewidth = 1)
		# ax2.set_ylabel("Ribosome Elongation Rate (aa/s)", fontsize=FONT_SIZE)
		# ax2.set_ylim([0,25])

		# # ax3.plot(doublingPerHour, stableRnaFraction, color = 'k', linewidth = 1)
		# # ax3.errorbar(doublingPerHour, stableRnaFraction, yerr=stableRnaFraction_error, color = "k", linewidth = 1)
		# # ax3.plot(doublingPerHourDB, stableRnaFractionDB, color = 'blue', linewidth = 1)
		# # ax3.errorbar(doublingPerHourDB, stableRnaFractionDB, yerr = stableRnaFractionDB_error, color = "blue", linewidth=1)
		# # ax3.set_ylabel("Synthesis rate Stable RNA / Total RNA")
		# # ax3.set_ylim([0,1])

		# # ax4.plot(doublingPerHour,mRnaFractionCalc, color='k', linewidth=1)
		# # ax4.errorbar(doublingPerHour, mRnaFractionCalc, yerr=mRnaFractionCalc_error, color = "k", linewidth = 1)
		# # ax4.plot(doublingPerHourDB,mRnaFractionDB, color='blue', linewidth=1)
		# # ax4.errorbar(doublingPerHourDB,mRnaFractionDB, yerr = mRnaFractionDB_error, color="blue", linewidth=1)
		# # ax4.set_ylabel("Synthesis rate mRNA / Total RNA")
		# # ax4.set_ylim([1,0])

		# ax5.plot(doublingPerHour,numOriginsAtInit, color='k', linewidth=1, label = 'Simulation')
		# ax5.errorbar(doublingPerHour, numOriginsAtInit, yerr=numOriginsAtInit_error, color = "k", linewidth = 1)
		# ax5.plot(doublingPerHourDB,numOriginsAtInitDB, color='blue', linewidth=1, label = 'Dennis and Bremer')
		# ax5.errorbar(doublingPerHourDB,numOriginsAtInitDB, yerr=numOriginsAtInitDB_error, linewidth=1, color="blue")
		# ax5.set_ylabel("Origins / Cell at Initiation", fontsize=FONT_SIZE)
		# ax5.set_ylim([0,5])
		# ax5.legend(loc=4,frameon=False, fontsize=FONT_SIZE)


		# whitePadSparklineAxis(ax0, False)
		# #whitePadSparklineAxis(ax1, False)
		# whitePadSparklineAxis(ax2, False)
		# #whitePadSparklineAxis(ax3, False)
		# #whitePadSparklineAxis(ax4)
		# whitePadSparklineAxis(ax5)

		# for a in [ax0, ax2, ax5]:
		# 	for tick in a.yaxis.get_major_ticks():
		# 		tick.label.set_fontsize(FONT_SIZE)
		# 	for tick in a.xaxis.get_major_ticks():
		# 		tick.label.set_fontsize(FONT_SIZE)

		# #ax4.set_xlabel("Doublings per Hour")
		# ax5.set_xlabel("Doublings per Hour", fontsize=FONT_SIZE)


		# #plt.subplots_adjust(wspace=0.25, hspace=0.25)

		# exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		# plt.close("all")


if __name__ == "__main__":
	Plot().cli()
