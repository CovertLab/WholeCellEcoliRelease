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

def main(seedOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

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
				#"proteinMass",
				#"tRnaMass",
				#"rRnaMass",
				#'mRnaMass',
				#"dnaMass"
				]

	cleanNames = [
				"Dry\nmass",
				#"Protein\nmass",
				#"tRNA\nmass",
				#"rRNA\nmass",
				#"mRNA\nmass",
				#"DNA\nmass"
				]

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx])

			f = plt.figure(figsize = (1.25, 0.8), frameon = False)
			#f = plt.figure(figsize = (4, 3), frameon = False)
			ax = f.add_axes([0, 0, 1, 1])
			ax.axis("off")

			ax.plot(time, massToPlot, linewidth = 2)
			#plt.plot(time / 60. / 60., massToPlot, linewidth = 2)
			ax.set_ylim([massToPlot.min(), massToPlot.max()])
			ax.set_xlim([time.min(), time.max()])
			#plt.ylim([massToPlot.min(), massToPlot.max()])
			#plt.xlim([time.min(), time.max()])
			print [massToPlot.min(), massToPlot.max()]

			from wholecell.analysis.analysis_tools import exportFigure
			exportFigure(plt, plotOutDir, "r01_{}_gen{}".format(massType, allDir.index(simDir)))
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
