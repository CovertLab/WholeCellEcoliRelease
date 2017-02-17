#!/usr/bin/env python

"""
Plot to compare cell properties across different growth conditions similar to Dennis and Bremer. 1996. Fig 2
Multiple generations of the same variant will be plotted together

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/7/16
"""

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from wholecell.utils.sparkline import whitePadSparklineAxis

NT_MW = 487.0
PROTEIN_MW = 110.0

FONT_SIZE=9

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):

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

	simDataFile = ap.get_variant_kb(all_cells[0])
	sim_data = cPickle.load(open(simDataFile, "rb"))
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

	#comparisonArray = np.zeros((4,3))

	#comparisonArray[0:] = rnaToProteinDB/rnaToProtein
	#comparisonArray[1:] = dnaToProteinDB/dnaToProtein
	#comparisonArray[2:] = elngRateDB/elngRate
	#comparisonArray[3:] = stableRnaFractionDB/stableRnaFraction
	#omparisonArray[4:]	=

	#comparisons

	# fig, ax = plt.subplots(3,2)#, figsize = (8.5, 11))

	# plt.subplots_adjust(wspace=0.4, hspace=0.4)

	# ax[0,0].scatter(rnaToProteinDB, rnaToProtein)
	# maxValue = np.max([max(rnaToProteinDB), max(rnaToProtein)])
	# ax[0,0].plot([0, maxValue],[0, maxValue])
	# ax[0,0].set_title("RNA to Protein (nuc/100 aa)", fontsize = FONT_SIZE)


	# ax[0,1].scatter(dnaToProteinDB, dnaToProtein)
	# maxValue = np.max([max(dnaToProteinDB), max(dnaToProtein)])
	# ax[0,1].plot([0, maxValue],[0, maxValue])
	# ax[0,1].set_title("DNA to Protein (chrom eq/10^9 aa)", fontsize = FONT_SIZE)


	# ax[1,0].scatter(elngRateDB, elngRate)
	# maxValue = np.max([max(elngRateDB), max(elngRate)])
	# ax[1,0].plot([0, maxValue],[0, maxValue])
	# ax[1,0].set_title("Ribosome Elongation Rate (aa/s)", fontsize = FONT_SIZE)

	# ax[1,1].scatter(stableRnaFractionDB, stableRnaFraction)
	# maxValue = np.max([max(stableRnaFractionDB), max(stableRnaFraction)])
	# ax[1,1].plot([0, maxValue],[0, maxValue])
	# ax[1,1].set_title("Rate Stable RNA to Rate Total RNA", fontsize = FONT_SIZE)

	# ax[2,0].scatter(mRnaFractionDB, mRnaFractionCalc)
	# maxValue = np.max([max(mRnaFractionDB), max(mRnaFractionCalc)])
	# ax[2,0].plot([0, maxValue],[0, maxValue])
	# ax[2,0].set_title("Rate mRNA to Rate Total RNA", fontsize = FONT_SIZE)

	# ax[2,1].scatter(numOriginsAtInitDB, numOriginsAtInit)
	# maxValue = np.max([max(numOriginsAtInitDB), max(numOriginsAtInit)])
	# ax[2,1].plot([0, maxValue],[0, maxValue])
	# ax[2,1].set_title("Origins per Cell at Initiation", fontsize = FONT_SIZE)

	# for a in ax.flatten():
	# 	whitePadSparklineAxis(a)
	# 	a.set_xlabel("Target value", fontsize = FONT_SIZE)
	# 	a.set_ylabel("Simulated value", fontsize = FONT_SIZE)

	fig, ax = plt.subplots(3,2, sharex=True, figsize = (10, 12))

	ax0 = ax[0,0]
	ax1 = ax[0,1]
	ax2 = ax[1,0]
	ax3 = ax[1,1]
	ax4 = ax[2,0]
	ax5 = ax[2,1]

	ax0.plot(doublingPerHour, rnaToProtein, color = 'k', linewidth = 2)
	ax0.errorbar(doublingPerHour, rnaToProtein, yerr=rnaToProtein_error, color = "k", linewidth = 2)
	ax0.plot(doublingPerHourDB, rnaToProteinDB, color = 'blue', linewidth = 2)
	ax0.errorbar(doublingPerHourDB, rnaToProteinDB, yerr = rnaToProteinDB_error, color = "blue", linewidth = 2)
	ax0.set_ylabel("RNA / Protein (nuc/100 aa)")
	ax0.set_ylim([0,20])

	ax1.plot(doublingPerHour, dnaToProtein, color = 'k', linewidth = 2)
	ax1.errorbar(doublingPerHour, dnaToProtein, yerr=dnaToProtein_error, color = "k", linewidth = 2)
	ax1.plot(doublingPerHourDB, dnaToProteinDB, color = 'blue', linewidth = 2)
	ax1.errorbar(doublingPerHourDB, dnaToProteinDB, yerr=dnaToProteinDB_error, color="blue", linewidth = 2)
	ax1.set_ylabel("DNA / Protein (chrom eq/" + r"$10^9$" + " aa)")
	ax1.set_ylim([0,4])

	ax2.plot(doublingPerHour, elngRate, color = 'k', linewidth = 2)
	ax2.errorbar(doublingPerHour, elngRate, yerr=elngRate_error, color = "k", linewidth = 2)
	ax2.plot(doublingPerHourDB, elngRateDB, color = 'blue', linewidth = 2)
	ax2.errorbar(doublingPerHourDB, elngRateDB, yerr = elngRateDB_error, color = "blue", linewidth = 2)
	ax2.set_ylabel("Ribosome Elongation Rate (aa/s)")
	ax2.set_ylim([0,25])

	ax3.plot(doublingPerHour, stableRnaFraction, color = 'k', linewidth = 2)
	ax3.errorbar(doublingPerHour, stableRnaFraction, yerr=stableRnaFraction_error, color = "k", linewidth = 2)
	ax3.plot(doublingPerHourDB, stableRnaFractionDB, color = 'blue', linewidth = 2)
	ax3.errorbar(doublingPerHourDB, stableRnaFractionDB, yerr = stableRnaFractionDB_error, color = "blue", linewidth=2)
	ax3.set_ylabel("Synthesis rate Stable RNA / Total RNA")
	ax3.set_ylim([0,1])

	ax4.plot(doublingPerHour,mRnaFractionCalc, color='k', linewidth=2) 
	ax4.errorbar(doublingPerHour, mRnaFractionCalc, yerr=mRnaFractionCalc_error, color = "k", linewidth = 2)
	ax4.plot(doublingPerHourDB,mRnaFractionDB, color='blue', linewidth=2)
	ax4.errorbar(doublingPerHourDB,mRnaFractionDB, yerr = mRnaFractionDB_error, color="blue", linewidth=2)
	ax4.set_ylabel("Synthesis rate mRNA / Total RNA")
	ax4.set_ylim([1,0])

	ax5.plot(doublingPerHour,numOriginsAtInit, color='k', linewidth=2, label = 'Simulation')
	ax5.errorbar(doublingPerHour, numOriginsAtInit, yerr=numOriginsAtInit_error, color = "k", linewidth = 2)
	ax5.plot(doublingPerHourDB,numOriginsAtInitDB, color='blue', linewidth=2, label = 'Dennis and Bremer')
	ax5.errorbar(doublingPerHourDB,numOriginsAtInitDB, yerr=numOriginsAtInitDB_error, linewidth=2, color="blue")
	ax5.set_ylabel("Origins / Cell at Initiation")
	ax5.set_ylim([0,5])
	ax5.legend(loc=4,frameon=False, fontsize=FONT_SIZE)


	whitePadSparklineAxis(ax0, False)
	whitePadSparklineAxis(ax1, False)
	whitePadSparklineAxis(ax2, False)
	whitePadSparklineAxis(ax3, False)
	whitePadSparklineAxis(ax4)
	whitePadSparklineAxis(ax5)

	ax4.set_xlabel("Doublings per Hour")
	ax5.set_xlabel("Doublings per Hour")


	plt.subplots_adjust(wspace=0.25, hspace=0.25)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")
	
	

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])