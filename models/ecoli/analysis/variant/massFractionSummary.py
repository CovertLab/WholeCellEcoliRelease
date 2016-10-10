#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import itertools

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from wholecell.analysis.plotting_tools import COLORS_LARGE

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):

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

	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	ap = AnalysisPaths(inputDir, variant_plot = True)
	all_cells = ap.get_cells()

	# Build a mapping from variant id to color
	idToColor = {}
	for idx, (cell_id, color) in enumerate(itertools.izip(all_cells, itertools.cycle(COLORS_LARGE))):
		idToColor[idx] = color

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fig, axesList = plt.subplots(len(massNames), sharex = True)

	currentMaxTime = 0
	for cellIdx, simDir in enumerate(all_cells):
		with open(os.path.join(simDir[:-32],'metadata','short_name')) as file:
			variant_name = [line for line in file][0]

		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for massIdx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massType)
			axesList[massIdx].plot(((time / 60.) / 60.), massToPlot, linewidth = 2, color=idToColor[cellIdx], label=variant_name)

			# set axes to size that shows all generations
			cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
			if cellCycleTime > currentMaxTime:
				currentMaxTime = cellCycleTime

			axesList[massIdx].set_xlim(0, currentMaxTime*ap.n_generation*1.1)
			axesList[massIdx].set_ylabel(cleanNames[massIdx] + " (fg)")

	for idx, axes in enumerate(axesList):
		axes.get_ylim()
		axes.set_yticks(list(axes.get_ylim()))

	axesList[0].set_title("Cell mass fractions")
	plt.legend(bbox_to_anchor=(.92, 5), loc=2, borderaxespad=0., prop={'size':6})
	axesList[len(massNames) - 1].set_xlabel("Time (hr)")
	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

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
