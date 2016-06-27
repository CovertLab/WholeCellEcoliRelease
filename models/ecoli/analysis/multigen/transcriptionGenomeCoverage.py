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

	mRnaIds = np.where(isMRna)[0]
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# Get number of mRNAs transcribed
	numGenesTranscribed = []
	mRnasTranscribedPerGen = []
	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIds]
		bulkMolecules.close()

		moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
		mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])
		numGenesTranscribed.append(sum(mRnasTranscribed))

		# mRnasNotTranscribed = np.logical_not(mRnasTranscribed)
		mRnasTranscribedNames = [mRnaNames[x] for x in np.arange(mRnaNames.shape[0]) if mRnasTranscribed[x]]
		mRnasTranscribedPerGen.append(mRnasTranscribedNames)

	# Plot
	numGens = len(numGenesTranscribed)
	ax = plt.subplot(1, 1, 1)
	ax.scatter(np.arange(numGens), numGenesTranscribed)	
	ax.set_xlabel("Generation", fontsize = 10)
	ax.set_ylabel("mRNAs transcribed", fontsize = 10)
	ax.set_title("Number of mRNAs transcribed\n(" + str(mRnaIds.shape[0]) + " total mRNA transcripts)", fontsize = 10)
	ax.set_xticks(np.arange(numGens))
	ax.set_ylim([0, 4353])

	# Plot names of transcribed mRNAs on the plot

	# for gen in np.arange(numGens):
	# 	text = str()
	# 	for mRna in mRnasTranscribedPerGen[gen]:
	# 		text += str(mRna) + "\n"
	# 	ax.text(gen, numGenesTranscribed[gen], text, fontsize = 6)


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


