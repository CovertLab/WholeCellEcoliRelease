#!/usr/bin/env python

from __future__ import division

import argparse
import os

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

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

AXIS_PADDING = 0.1

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	isMRna = kb.rnaData["isMRna"]

	rnaIds = kb.rnaData["id"][isMRna]

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("moleculeIDs")

	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]

	bulkMolecules.close()

	expectedCountsArbitrary = kb.rnaExpression['expression'][isMRna]

	expectedFrequency = expectedCountsArbitrary/expectedCountsArbitrary.sum()

	actualFrequency = 1.* rnaCountsBulk.T / rnaCountsBulk.sum(1)

	plt.figure(figsize = (8.5, 11))

	nColors = len(COLORS)
	nTimesteps = actualFrequency.shape[1]
	onesVector = np.ones(nTimesteps)

	sorting = np.argsort(expectedFrequency)

	expectedFrequency = expectedFrequency[sorting]
	actualFrequency = actualFrequency[sorting, :]

	maxValue = max(expectedFrequency.max(), actualFrequency.max()) * (1 + AXIS_PADDING)

	plt.plot([0, maxValue], [0, maxValue], 'k--', linewidth = 0.2)

	# Plot style inspired by the prettyplotlib package (see:
	# http://blog.olgabotvinnik.com/post/58941062205/prettyplotlib-painlessly-create-beautiful-matplotlib)

	colors = np.array(
		(COLORS*(expectedFrequency.size // len(COLORS) + 1))[:expectedFrequency.size]
		).repeat(actualFrequency.shape[1], 0)

	plt.scatter(
		expectedFrequency.repeat(actualFrequency.shape[1]),
		actualFrequency.reshape(-1),
		c = colors,
		edgecolors = "none",
		rasterized = True,
		alpha = 0.1
		)

	plt.scatter(
		expectedFrequency, actualFrequency[:, -1],
		c = "none",
		edgecolors = "k",
		rasterized = True
		)

	plt.xlabel("Expected mRNA frequency")
	plt.ylabel("Actual mRNA frequency (all time steps)")
	plt.title("Simulation gene expression")

	plt.axis("scaled")
	plt.axis([0, maxValue]*2)

	# plt.show()

	# plotOutFileName = plotOutFileName[:plotOutFileName.index(".")] + ".png" # hack to save this as raster

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


if __name__ == "__main__":
	# defaultKBFile = os.path.join(
	# 		wholecell.utils.constants.SERIALIZED_KB_DIR,
	# 		wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
	# 		)

	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	# parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])

	print "Disabled because it's slow"
