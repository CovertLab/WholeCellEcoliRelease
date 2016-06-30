#!/usr/bin/env python
"""
Plots fraction of mRNAs transcribed (out of all genes to be transcribed) for all generations.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2016
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

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	# Get IDs of mRNAs
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	degRate = sim_data.process.transcription.rnaData["degRate"]
	basalExpression = sim_data.process.transcription.rnaExpression["basal"]
	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]
	mRnaIds = np.where(isMRna)[0]

	mRnaBasalExpression = np.array([basalExpression[x] for x in mRnaIds])
	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIds])
	mRnaDegRate = np.array([degRate[x] for x in mRnaIds])
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# Sort in order of decreasing basal expression
	descendingOrderIndexing = np.argsort(mRnaBasalExpression)[::-1]
	mRnaBasalExpressionSorted = mRnaBasalExpression[descendingOrderIndexing]
	mRnaSynthProbSorted = mRnaSynthProb[descendingOrderIndexing]
	mRnaDegRateSortedUnits = mRnaDegRate[descendingOrderIndexing]
	mRnaNamesSorted = mRnaNames[descendingOrderIndexing]

	# Remove units on degradation rates
	mRnaDegRateSorted = np.array([x.asNumber(1 / units.s) for x in mRnaDegRateSortedUnits])

	mRnaDataSorted = {	0: {"name": "Basal expression", "data": mRnaBasalExpressionSorted, "zoomMax": 0.000005}, 
						1: {"name": "Synthesis probability", "data": mRnaSynthProbSorted, "zoomMax": 0.0001},
						2: {"name": "Degradation rate", "data": mRnaDegRateSorted, "zoomMax": 0.005},
					}

	# Get proteins of interest
	dcuRId = np.where([x == "G7826_RNA[c]" for x in mRnaNamesSorted])[0]
	baeRId = np.where([x == "EG11618_RNA[c]" for x in mRnaNamesSorted])[0]
	narLId = np.where([x == "EG10643_RNA[c]" for x in mRnaNamesSorted])[0]

	proteinsOfInterest = np.array([dcuRId, baeRId, narLId])
	proteinsOfInterestNames = np.array(["dcuR", "baeR", "narL"])
	
	# Get number of mRNAs transcribed
	transcribedBool = []

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNamesSorted])
		moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]
		bulkMolecules.close()

		moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
		mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])

		mRnasTranscribedNames = [mRnaNames[x] for x in np.arange(mRnaNames.shape[0]) if mRnasTranscribed[x]]
		transcribedBool.append(mRnasTranscribed)

	transcribedBool = np.array(transcribedBool)

	# Plot
	numGens = allDir.shape[0]
	numDataSorted = len(mRnaDataSorted)
	numMRnas = mRnaNamesSorted.shape[0]
	fig = plt.figure(figsize = (15, 18))
	border = 0

	rows = numGens + 2*numDataSorted + 1
	cols = 1
	xvals = np.arange(numMRnas)
	intersection = np.ones((1, numMRnas), dtype = bool)


	# Plot sorted mRNA data as top 6 subplots
	for idx in np.arange(2*numDataSorted):
		ax = plt.subplot(rows, cols, idx + 1)
		ax.vlines(xvals, [0], mRnaDataSorted[idx / 2]["data"], colors = "0.15")
		ax.set_xlim([-border, numMRnas + border])
		ax.set_ylim([np.min(mRnaDataSorted[idx / 2]["data"]), np.max(mRnaDataSorted[idx / 2]["data"])])
		ax.tick_params(which = "both", direction = "out", right = "off", top = "off")
		ax.spines["left"].set_visible(False)
		ax.spines["right"].set_visible(False)

		heighOffset = np.max(mRnaDataSorted[idx / 2]["data"]) / 5.
		for proteinIdx, protein in enumerate(proteinsOfInterest):
			ax.vlines(protein, [0], np.max(mRnaDataSorted[idx / 2]["data"]), color = "orange")
			ax.annotate(proteinsOfInterestNames[proteinIdx], xy = (protein, np.max(mRnaDataSorted[idx / 2]["data"]) - (proteinIdx + 1) * heighOffset))

		if idx %2 == 0:
			ax.set_title(mRnaDataSorted[idx / 2]["name"])
			ax.set_yticks([np.min(mRnaDataSorted[idx / 2]["data"]), np.max(mRnaDataSorted[idx / 2]["data"])])

		else:
			ax.set_ylim([np.min(mRnaDataSorted[idx / 2]["data"]), mRnaDataSorted[idx / 2]["zoomMax"]])
			ax.set_yticks([np.min(mRnaDataSorted[idx / 2]["data"]), mRnaDataSorted[idx / 2]["zoomMax"]])


	# Plot binary plot for each generation
	for idx, subplotIdx in enumerate(np.arange(7, 7 + numGens)):
		ax = plt.subplot(rows, cols, subplotIdx)
		ax.vlines(xvals, [0], transcribedBool[idx], color = "0.15")
		ax.set_title("Generation %s" % idx, fontsize = 14)
		ax.tick_params(which = "both", direction = "out", right = "off", top = "off")
		ax.set_yticks([])
		ax.set_xlim([-border, numMRnas + border])
		ax.spines["left"].set_visible(False)
		ax.spines["right"].set_visible(False)

		for protein in proteinsOfInterest:
			if transcribedBool[idx][protein]:
				ax.plot(np.ones(2)*protein, [0, 1], linestyle = "--", color = "g")
			else:
				ax.plot(np.ones(2)*protein, [0, 1], linestyle = "--", color = "r")

		intersection = np.logical_and(intersection, transcribedBool[idx])


	# Plot intersection binary plot
	ax = plt.subplot(rows, cols, rows)
	highlightWhitespaceBoundary = 2500
	ax.vlines(xvals[:highlightWhitespaceBoundary], [0], intersection[0][:highlightWhitespaceBoundary])
	ax.vlines(xvals, [0], np.logical_not(intersection), color = "#3399ff")
	ax.vlines(xvals[highlightWhitespaceBoundary:], [0], intersection[0][highlightWhitespaceBoundary:])
	
	ax.set_title("Intersection of all generations", fontsize = 14)
	ax.set_xlabel("mRNA transcripts\n(in order of decreasing expected basal expression)", fontsize = 10)
	ax.set_yticks([])
	ax.set_xlim([-border, numMRnas + border])
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)


	plt.subplots_adjust(hspace = 1, wspace = 0)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	plt.close("all")


	# # Plot polar coordinate representation of the most differing set of mRNAs across generations
	# mostDifferentlyTranscribed = []
	# for mRnaIdx in np.arange(numMRnas):
	# 	if np.sum(transcribedBool[:, mRnaIdx]) not in [0, numGens]:
	# 		mostDifferentlyTranscribed.append(mRnaIdx)

	# mostDifferentlyTranscribed = np.array(mostDifferentlyTranscribed)
	# numMRnasPolar = mostDifferentlyTranscribed.shape[0]

	# fig = plt.figure(figsize = (12, 12))
	# plotOutFileName = "transcriptionGenomeCoveragePolar"

	# rGen = 1
	# rSpace = 0.25
	# ax.set_ylim([0, ((rSpace + rSpace)* numGens)])
	# ax = plt.subplot(1, 1, 1, polar = True)
	# colors = ["0.2", "0.4", "0.6", "0.8"]

	# for idx in np.arange(numGens - 1, -1, -1):
	# 	radius = np.ones(numMRnasPolar) * (idx + 1) * (rGen + rSpace) * transcribedBool[idx, mostDifferentlyTranscribed]
	# 	ax.plot(np.linspace(0, 2*np.pi, numMRnasPolar), radius, colors[idx])
	# 	ax.fill_between(np.linspace(0, 2*np.pi, numMRnasPolar), radius, 0, color = colors[idx])

	# 	# white space
	# 	# rWhiteSpace = np.ones(numMRnasPolar) * idx * (rGen + rSpace)
	# 	# ax.fill_between(np.linspace(0, 2*np.pi, numMRnasPolar), rWhiteSpace + rSpace, 0, color = "white")

	# radius = numGens * (rGen + rSpace)
	# thetaVals = np.linspace(0, 2*np.pi, numMRnasPolar)
	# for thetaIndex, idx in enumerate(mostDifferentlyTranscribed):
	# 	ax.annotate(str(mRnaNamesSorted[mostDifferentlyTranscribed[thetaIndex]]),
	# 				xy = (thetaVals[thetaIndex], radius),
	# 				)
	# ax.grid(False)
	# ax.set_xticks([])
	# ax.set_yticks([])

	# from wholecell.analysis.analysis_tools import exportFigure
	# exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	# plt.close("all")

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


