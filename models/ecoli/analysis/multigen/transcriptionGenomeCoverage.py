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
	basalExpression = sim_data.process.transcription.rnaExpression["basal"]
	# basalExpression = sim_data.process.transcription.rnaSynthProb["basal"]



	mRnaIds = np.where(isMRna)[0]
	mRnaBasalExpression = np.array([basalExpression[x] for x in mRnaIds])
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# Sort in order of decreasing basal expression
	descendingOrderIndexing = np.argsort(mRnaBasalExpression)[::-1]
	mRnaBasalExpressionSorted = mRnaBasalExpression[descendingOrderIndexing]
	mRnaNamesSorted = mRnaNames[descendingOrderIndexing]

	# Get number of mRNAs transcribed
	numGenesTranscribed = []
	mRnasTranscribedPerGen = []
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
		numGenesTranscribed.append(np.sum(mRnasTranscribed))

		mRnasTranscribedNames = [mRnaNames[x] for x in np.arange(mRnaNames.shape[0]) if mRnasTranscribed[x]]
		mRnasTranscribedPerGen.append(mRnasTranscribedNames)

		transcribedBool.append(mRnasTranscribed)

	transcribedBool = np.array(transcribedBool)
	numGenesTranscribed = np.array(numGenesTranscribed)

	# Plot
	numGens = numGenesTranscribed.shape[0]
	numMRnas = mRnaNamesSorted.shape[0]
	fig = plt.figure(figsize = (20, 10))

	rows = numGens + 1
	cols = 1
	xvals = np.arange(numMRnas)
	intersection = np.ones((1, numMRnas), dtype = bool)

	for idx in np.arange(numGens):
		ax = plt.subplot(rows, cols, idx + 1)
		ax.vlines(xvals, [0], transcribedBool[idx], color = "0.25")
		ax.set_title("Generation %s" % idx, fontsize = 14)
		ax.set_yticks([])
		ax.set_xlim([0, numMRnas])

		ax.vlines(xvals[:2600], [0], np.logical_not(transcribedBool[idx][:2600]), color = "#3399ff")

		intersection = np.logical_and(intersection, transcribedBool[idx])

	ax = plt.subplot(rows, cols, rows)
	ax.set_title("Intersection of all generations", fontsize = 14)
	ax.vlines(xvals, [0], intersection[0])
	# ax.vlines(xvals[:2500], [0], np.logical_not(intersection[idx][:2500]), color = "#3399ff")
	ax.set_xlabel("mRNA transcripts\nin order of decreasing expected basal expression", fontsize = 14)
	ax.set_yticks([])
	ax.set_xlim([0, numMRnas])

	plt.subplots_adjust(hspace = 0.5, wspace = 0)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	plt.close("all")


	# Plot polar representation of the most differing set of mRNAs
	mostDifferentlyTranscribed = []
	for mRnaIdx in np.arange(numMRnas):
		if np.sum(transcribedBool[:, mRnaIdx]) not in [0, 3]:
			mostDifferentlyTranscribed.append(mRnaIdx)

	mostDifferentlyTranscribed = np.array(mostDifferentlyTranscribed)
	numMRnasPolar = mostDifferentlyTranscribed.shape[0]

	fig = plt.figure(figsize = (12, 12))
	plotOutFileName = "transcriptionGenomeCoveragePolar"

	rGen = 1
	rSpace = 0.25
	ax.set_ylim([0, ((rSpace + rSpace)* numGens)])
	ax = plt.subplot(1, 1, 1, polar = True)

	for idx in np.arange(numGens - 1, -1, -1):
		radius = np.ones(numMRnasPolar) * (idx + 1) * (rGen + rSpace) * transcribedBool[idx, mostDifferentlyTranscribed]
		ax.plot(np.linspace(0, 2*np.pi, numMRnasPolar), radius)
		ax.fill_between(np.linspace(0, 2*np.pi, numMRnasPolar), radius, 0)

	ax.plot(np.linspace(0, 2*np.pi, 2), [True, False])

	ax.grid(False)
	ax.set_xticks([])
	ax.set_yticks([])

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


