#!/usr/bin/env python
"""
Plots transcript frequency (ie. frequency of observing at least 
one copy of transcript) at the 4th generation across 32 seeds.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/27/2017
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

N_SEEDS = 32

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	return
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	allDir = ap.get_cells(generation = [3])

	if len(allDir) <= 1:
		print "Skipping -- transcriptFrequency only runs for multiple seeds"
		return

	sim_data = cPickle.load(open(simDataFile, "rb"))

	# Get mRNA data
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	mRnaIndexes = np.where(isMRna)[0]
	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]

	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIndexes])
	mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

	# Get mRNA indices in bulk molecules
	simOutDir = os.path.join(allDir[0], "simOut")
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	mRnaIndexes_bulk = [moleculeIds.index(x) for x in mRnaIds]
	bulkMolecules.close()

	# Get frequency data
	hadTranscribed = None
	simulatedSynthProb = None
	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes_bulk]
		bulkMolecules.close()
		moleculeCounts_sumOverTime = moleculeCounts.sum(axis = 0)
		mRnasTranscribed = np.array([x != 0 for x in moleculeCounts_sumOverTime])

		rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		simulatedSynthProb_ = rnaSynthProb.readColumn("rnaSynthProb")
		rnaSynthProb.close()

		if hadTranscribed == None:
			hadTranscribed = mRnasTranscribed
			simulatedSynthProb = np.mean(simulatedSynthProb_, axis = 0)[isMRna]
			continue
		hadTranscribed = np.vstack((hadTranscribed, mRnasTranscribed))
		simulatedSynthProb = np.vstack((simulatedSynthProb, np.mean(simulatedSynthProb_, axis = 0)[isMRna]))

		if hadTranscribed.shape[0] == N_SEEDS:
			break

	hadTranscribedFrequency = np.mean(hadTranscribed, axis = 0)
	simulatedSynthProbAvg = np.mean(simulatedSynthProb, axis = 0)
	indexOrder = np.argsort(simulatedSynthProbAvg)

	# Plot
	fig = plt.figure(figsize = (16, 8))
	scatterAxis = plt.subplot2grid((2, 4), (0, 0), colspan = 3, rowspan = 2)
	histAxis = plt.subplot2grid((2, 4), (0, 3), colspan = 1, rowspan = 2, sharey = scatterAxis)

	colors = np.repeat("g", len(hadTranscribedFrequency))
	colors[np.where(hadTranscribedFrequency[indexOrder] == 1)] = "b"
	colors[np.where(hadTranscribedFrequency[indexOrder] == 0)] = "r"

	scatterAxis.scatter(np.arange(len(hadTranscribedFrequency)), hadTranscribedFrequency[indexOrder], marker = "o", facecolors = colors, edgecolors = "none", s = 20)
	scatterAxis.set_xlim([0, len(hadTranscribedFrequency)])
	scatterAxis.set_ylim([0, 1])
	whitePadSparklineAxis(scatterAxis)

	N, bins, patches = histAxis.hist(hadTranscribedFrequency, bins = N_SEEDS + 1, orientation = 'horizontal')
	for i in xrange(1, len(patches) - 1):
		plt.setp(patches[i], facecolor = "none", edgecolor = "g")
	plt.setp(patches[0], facecolor = "none", edgecolor = "r")
	plt.setp(patches[-1], facecolor = "none", edgecolor = "b")
	histAxis.set_xscale("log")
	whitePadSparklineAxis(histAxis)
	histAxis.xaxis.tick_bottom()
	histXmin, histXmax = histAxis.get_xlim()

	scatterAxis.set_ylim([-.01, 1.01])
	plt.subplots_adjust(wspace = 0.6, hspace = 0.4, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close()

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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])