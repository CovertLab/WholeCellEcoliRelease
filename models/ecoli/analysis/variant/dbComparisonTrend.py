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

NT_MW = 487.0
PROTEIN_MW = 110.0

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):

	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	ap = AnalysisPaths(inputDir, variant_plot = True)
	all_cells = ap.get_cells()

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
	#import ipdb
	#ipdb.set_trace()


	#import all Dennis and Bremer Data:
	#dt order 0.6, 1.5, 2.5
	rnaToProteinDB = [6.615, 10.962, 15.6]
	dnaToProteinDB = [2.82, 1.79, 1.52]
	elngRateDB = [12, 18, 21]
	stableRnaFractionDB = [.41,.68,.85]
	mRnaFractionDB = [.59, .32, .15]
	stableRnaSynthRateDB = [0.53, 2.23, 5.3]
	numOriginsAtInitDB = [1,2,4]
	doublingPerHourDB = [0.6, 1.5, 2.5]

	#comparisonArray = np.zeros((4,3))

	#comparisonArray[0:] = rnaToProteinDB/rnaToProtein
	#comparisonArray[1:] = dnaToProteinDB/dnaToProtein
	#comparisonArray[2:] = elngRateDB/elngRate
	#comparisonArray[3:] = stableRnaFractionDB/stableRnaFraction
	#omparisonArray[4:]	=

	#comparisons


	fig, ax = plt.subplots(5, sharex=True, figsize = (8.5, 11))

	ax[0].plot(doublingPerHour, rnaToProtein, color = 'k', linewidth = 2)
	ax[0].plot(doublingPerHourDB, rnaToProteinDB, color = 'r', linewidth = 2)
	ax[0].set_ylabel("RNA to Protein\n(nuc/100 aa)")
	ax[0].set_ylim([0,20])

	ax[1].plot(doublingPerHour, dnaToProtein, color = 'k', linewidth = 2)
	ax[1].plot(doublingPerHourDB, dnaToProteinDB, color = 'r', linewidth = 2)
	ax[1].set_ylabel("DNA to Protein\n(chrom eq/10^9 aa)")
	ax[1].set_ylim([0,4])

	ax[2].plot(doublingPerHour, elngRate, color = 'k', linewidth = 2)
	ax[2].plot(doublingPerHourDB, elngRateDB, color = 'r', linewidth = 2)
	ax[2].set_ylabel("Ribosome Elongation\nRate (aa/s)")
	ax[2].set_ylim([0,25])

	ax[3].plot(doublingPerHour, stableRnaFraction, color = 'k', linewidth = 2)
	ax[3].plot(doublingPerHourDB, stableRnaFractionDB, color = 'r', linewidth = 2)
	ax[3].set_ylabel("Rate Stable RNA to\nRate Total RNA")
	ax[3].set_ylim([0,1])

	b = ax[3].twinx()
	b.plot(doublingPerHour,mRnaFractionCalc, color='k', linewidth=2)
	b.plot(doublingPerHourDB,mRnaFractionDB, color='r', linewidth=2)
	b.set_ylabel("Rate mRNA to\nRate Total RNA")
	b.set_ylim([1,0])

	ax[4].plot(doublingPerHour,numOriginsAtInit, color='k', linewidth=2, label = 'Simulation')
	ax[4].plot(doublingPerHourDB,numOriginsAtInitDB, color='r', linewidth=2, label = 'Dennis and Bremer')
	ax[4].set_ylabel("Origins per Cell\nat Initiation")
	ax[4].set_xlabel("Doublings per Hour")
	ax[4].set_ylim([0,5])
	ax[4].legend(loc=4,frameon=False)


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
