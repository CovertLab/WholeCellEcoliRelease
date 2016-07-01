#!/usr/bin/env python
"""
Plots fraction of mRNAs transcribed (out of all genes to be transcribed) for all seeds.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/29/2016
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

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get IDs of mRNAs
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	basalExpression = sim_data.process.transcription.rnaExpression["basal"]
	mRnaIds = np.where(isMRna)[0]

	mRnaBasalExpression = np.array([basalExpression[x] for x in mRnaIds])
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# Sort in order of decreasing basal expression
	descendingOrderIndexing = np.argsort(mRnaBasalExpression)[::-1]
	mRnaBasalExpressionSorted = mRnaBasalExpression[descendingOrderIndexing]
	mRnaNamesSorted = mRnaNames[descendingOrderIndexing]

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells()

	# Get number of mRNAs transcribed
	transcribedFreq = []
	for simDir in all_cells:
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNamesSorted])
		moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]
		bulkMolecules.close()

		moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
		mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])

		transcribedFreq.append(mRnasTranscribed)

	transcribedFreq = np.array(transcribedFreq)
	transcribedFreqSumOverSeeds = transcribedFreq.sum(axis = 0)

	# Get proteins of interest
	dcuRId = np.where([x == "G7826_RNA[c]" for x in mRnaNamesSorted])[0]
	baeRId = np.where([x == "EG11618_RNA[c]" for x in mRnaNamesSorted])[0]
	narLId = np.where([x == "EG10643_RNA[c]" for x in mRnaNamesSorted])[0]

	proteinsOfInterest = np.array([dcuRId, baeRId, narLId])
	proteinsOfInterestNames = np.array(["dcuR", "baeR", "narL"])

	# Plot
	numMRnas = mRnaNamesSorted.shape[0]
	numSeeds = all_cells.shape[0]
	fig = plt.figure(figsize = (14, 10))

	ax = plt.subplot(1, 1, 1)
	ax.scatter(np.arange(numMRnas), transcribedFreqSumOverSeeds / float(numSeeds), facecolors = "none", edgecolors = "b")

	heightOffset = 0.01
	for proteinIdx, protein in enumerate(proteinsOfInterest):
		ax.scatter(protein, transcribedFreqSumOverSeeds[protein] / float(numSeeds), facecolors = "orange", edgecolors = "none")
		ax.annotate(proteinsOfInterestNames[proteinIdx], xy = (protein, transcribedFreqSumOverSeeds[protein] / float(numSeeds) + heightOffset))

	ax.set_title("Frequency of producing at least 1 transcript\n(n = %s seeds)" % numSeeds, fontsize = 12)
	ax.set_xlabel("mRNA transcripts\n(in order of decreasing expected basal expression)", fontsize = 10)
	ax.set_xlim([0, numMRnas])
	ax.set_ylim([-.05, 1.05])
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)

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