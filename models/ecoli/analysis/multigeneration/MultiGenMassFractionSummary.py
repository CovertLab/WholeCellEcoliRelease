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

	# Get all cells
	allDir = ap.getAll()

	massNames = ["dryMass","proteinMass","tRnaMass","rRnaMass",'mRnaMass',"dnaMass"]

	plt.figure(figsize = (8.5, 11))

	for idx, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		#initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for idx2, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx2])
			# massToPlot = massToPlot / massToPlot[0]

			plt.subplot(6,1,idx2+1)
			plt.plot(time / 60., massToPlot, linewidth = 2)
			plt.ylabel(massNames[idx2])
			if idx2 == 0:
				plt.title("Cell mass fractions (fg)")

	plt.xlabel("Time (min)")

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
