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

def main(seedOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir)

	# TODO: Declutter Y-axis

	# Get all cells
	allDir = ap.getAll()

	massNames = [
				"dryMass",
				"proteinMass",
				#"tRnaMass",
				"rRnaMass",
				'mRnaMass',
				"dnaMass"
				]

	cleanNames = [
				"Dry\nmass",
				"Protein\nmass",
				#"tRNA\nmass",
				"rRNA\nmass",
				"mRNA\nmass",
				"DNA\nmass"
				]

	#plt.figure(figsize = (8.5, 11))
	fig, axesList = plt.subplots(len(massNames), sharex = True)

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")
		#initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx])
			# massToPlot = massToPlot / massToPlot[0]

			axesList[idx].plot(time / 60. / 60., massToPlot, linewidth = 2)
			axesList[idx].set_xlim(0, 5.1)
			axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

	for axes in axesList:
		axes.get_ylim()
		axes.set_yticks(list(axes.get_ylim()))

	axesList[0].set_title("Cell mass fractions")
	axesList[len(massNames) - 1].set_xlabel("Time (hr)")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
