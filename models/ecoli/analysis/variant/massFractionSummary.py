#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[228,26,28],
	[55,126,184],
	[77,175,74],
	[152,78,163],
	[255,127,0],
	[255,255,51],
	[166,86,40],
	[247,129,191]
	]

COLORS = [
	[colorValue/255. for colorValue in color]
	for color in COLORS_256
	]

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

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fig, axesList = plt.subplots(len(massNames), sharex = True)

	currentMaxTime = 0
	for simDir in all_cells:
		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massType)
			axesList[idx].plot(((time / 60.) / 60.), massToPlot, linewidth = 2)

			# set axes to size that shows all generations
			cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
			if cellCycleTime > currentMaxTime:
				currentMaxTime = cellCycleTime

			axesList[idx].set_xlim(0, currentMaxTime*int(metadata["total_gens"])*1.1)
			axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

	for axes in axesList:
		axes.get_ylim()
		axes.set_yticks(list(axes.get_ylim()))
	
	axesList[0].set_title("Cell mass fractions")
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
