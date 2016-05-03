#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# TODO: Declutter Y-axis

	# Get first cell from each generation
	firstCellLineage = []
	for gen_idx in range(ap.n_generation):
		firstCellLineage.append(ap.get_cells(generation = [gen_idx])[0])

	massNames = [
				#"dryMass",
				"proteinMass",
				"tRnaMass",
				"rRnaMass",
				'mRnaMass',
				"dnaMass"
				]

	cleanNames = [
				#"Dry\nmass",
				"Protein\nmass frac.",
				"tRNA\nmass frac.",
				"rRNA\nmass frac.",
				"mRNA\nmass frac.",
				"DNA\nmass frac."
				]

	fig, axesList = plt.subplots(len(massNames), sharex = True)

	for simDir in firstCellLineage:
		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		massData = np.zeros((len(massNames),time.size))

		for idx, massType in enumerate(massNames):
			massData[idx,:] = mass.readColumn(massNames[idx])

		massData = massData / massData.sum(axis = 0)

		for idx, massType in enumerate(massNames):
			axesList[idx].plot(time / 60, massData[idx,:])
			axesList[idx].set_ylabel(cleanNames[idx])

	for axes in axesList:
		axes.set_yticks(list(axes.get_ylim()))

	axesList[-1].set_xlabel('Time (min)')

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
