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

# TODO: account for complexation

CMAP_COLORS_255 = [
	[103,0,31],
	[178,24,43],
	[214,96,77],
	[244,165,130],
	[253,219,199],
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	rnaIds = kb.process.transcription.rnaData["id"][kb.relation.rnaIndexToMonomerMapping]

	proteinIds = kb.process.translation.monomerData["id"]

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]

	proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

	proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime	

	bulkMolecules.close()

	time /= 60.

	# TODO: plots for all molecules (needs to be fast)
	# TODO: account for complex composition, other modifications

	index = 611
	name = "DNA polymerase"

	plt.figure(figsize = (8.5, 11))

	# TODO: functionalize this as a filled step plot

	stepTime = np.roll(np.repeat(time, 2), +1)[1:]

	stepRna = np.repeat(rnaCountsBulk[:, index], 2)[1:]
	stepProtein = np.repeat(proteinCountsBulk[:, index], 2)[1:]

	plt.subplot(2, 1, 1)

	polygonRna = np.vstack([
		np.hstack([stepTime[0], stepTime, stepTime[-1]]),
		np.hstack([0, stepRna, 0])
		]).T

	plt.gca().add_patch(plt.Polygon(
		polygonRna,
		facecolor = CMAP_COLORS[3],
		edgecolor = "none"
		))

	lineRna = np.vstack([stepTime, stepRna]).T

	plt.gca().add_patch(plt.Polygon(
		lineRna,
		facecolor = "none",
		edgecolor = CMAP_COLORS[1],
		closed = False,
		linewidth = 2
		))

	plt.axis("auto")

	plt.xlabel("Time (min)")
	plt.ylabel("mRNA copies (number)")

	plt.subplot(2, 1, 2)

	polygonProtein = np.vstack([
		np.hstack([stepTime[0], stepTime, stepTime[-1]]),
		np.hstack([0, stepProtein, 0])
		]).T

	plt.gca().add_patch(plt.Polygon(
		polygonProtein,
		facecolor = CMAP_COLORS[-4],
		edgecolor = "none"
		))

	lineProtein = np.vstack([stepTime, stepProtein]).T

	plt.gca().add_patch(plt.Polygon(
		lineProtein,
		facecolor = "none",
		edgecolor = CMAP_COLORS[-2],
		closed = False,
		linewidth = 2
		))

	plt.axis("auto")

	plt.xlabel("Time (min)")
	plt.ylabel("Protein copies (number)")

	plt.suptitle(name)

	# plt.show()

	# import ipdb; ipdb.set_trace()

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
