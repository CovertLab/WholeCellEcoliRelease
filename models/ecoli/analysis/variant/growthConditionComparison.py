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

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)
	
	rnaToProteinDict = {}
	dnaToProteinDict = {}
	elngRateDict = {}
	stableRnaFractionDict = {}
	doublingPerHourDict = {}

	simDataFile = ap.get_variant_kb(all_cells[0])
	sim_data = cPickle.load(open(simDataFile, "rb"))
	nAvogadro = sim_data.constants.nAvogadro.asNumber()
	chromMass = (sim_data.getter.getMass(['CHROM_FULL[c]'])[0] / sim_data.constants.nAvogadro).asNumber() 
	
	for simDir in all_cells:
		simOutDir = os.path.join(simDir, "simOut")
		variant = int(simDir[simDir.rfind('generation_')-14:simDir.rfind('generation_')-8])

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
		stableRnaSynth = np.sum(rnaSynth[:,isTRna], axis=1) + np.sum(rnaSynth[:,isRRna], axis=1)
		totalRnaSynth = np.sum(rnaSynth, axis=1).astype(float)
		rnaFraction = stableRnaSynth / totalRnaSynth

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		uniqueMoleculeCounts.close()

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

		if variant in rnaToProteinDict.keys():
			rnaToProteinDict[variant] = np.append(rnaToProteinDict[variant], rnaNT / (proteinAA / 100))
			dnaToProteinDict[variant] = np.append(dnaToProteinDict[variant], chromEquivalents / (proteinAA / 10**9))
			elngRateDict[variant] = np.append(elngRateDict[variant], (actualElongations / activeRibosome / timeStepSec)[3:])
			stableRnaFractionDict[variant] = np.append(stableRnaFractionDict[variant], np.asarray(rnaFraction)[~np.isnan(rnaFraction)])
			doublingPerHourDict[variant] = np.append(doublingPerHourDict[variant], 60 / doublingTime)
		else:
			rnaToProteinDict[variant] = rnaNT / (proteinAA / 100)
			dnaToProteinDict[variant] = chromEquivalents / (proteinAA / 10**9)
			elngRateDict[variant] = (actualElongations / activeRibosome / timeStepSec)[3:]
			stableRnaFractionDict[variant] = np.asarray(rnaFraction)[~np.isnan(rnaFraction)]
			doublingPerHourDict[variant] = 60 / doublingTime

	rnaToProtein = []
	dnaToProtein = []
	elngRate = []
	stableRnaFraction = []
	doublingPerHour = []

	for key in rnaToProteinDict.keys():
		rnaToProtein += [rnaToProteinDict[key]]
		dnaToProtein += [dnaToProteinDict[key]]
		elngRate += [elngRateDict[key]]
		stableRnaFraction += [stableRnaFractionDict[key]]
		doublingPerHour += [np.mean(doublingPerHourDict[key])]

	plt.figure(figsize = (8.5, 11))

	sp = plt.subplot(4,1,1)
	sp.violinplot(rnaToProtein, positions=doublingPerHour, showmeans=True)
	sp.set_ylabel("RNA to Protein\n(nuc/100 aa)")

	sp = plt.subplot(4,1,2)
	sp.violinplot(dnaToProtein, positions=doublingPerHour, showmeans=True)
	sp.set_ylabel("DNA to Protein\n(chrom eq/10^9 aa)")

	sp = plt.subplot(4,1,3)
	sp.violinplot(elngRate, positions=doublingPerHour, showmeans=True)
	sp.set_ylabel("Ribosome Elongation\nRate (aa/s)")

	sp = plt.subplot(4,1,4)
	sp.violinplot(stableRnaFraction, positions=doublingPerHour, showmeans=True)
	sp.set_ylabel("Rate Stable RNA to\nRate Total RNA")
	sp.set_xlabel("Doublings per Hour")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

	raise Exception()	

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
