#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

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

	kb = cPickle.load(open(kbFile, "rb"))

	rnaIds = kb.process.transcription.rnaData["id"][kb.relation.rnaIndexToMonomerMapping]

	proteinIds = kb.process.translation.monomerData["id"]

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime	

	bulkMolecules.close()

	time /= 60.

	index = 611
	name = "DNA polymerase"

	stepTime = np.roll(np.repeat(time, 2), +1)[1:]

	stepRna = np.repeat(rnaCountsBulk[:, index], 2)[1:]

	lineRna = np.vstack([stepTime, stepRna]).T




	f = plt.figure(figsize = (1.25, 0.78), frameon = False)
	ax = f.add_axes([0, 0, 1, 1])
	ax.axis("off")

	ax.add_patch(plt.Polygon(
		lineRna,
		facecolor = "none",
		edgecolor = CMAP_COLORS[-2],
		closed = False,
		linewidth = 2
		))
	ax.set_ylim([lineRna[:, 1].min(), lineRna[:, 1].max()])
	ax.set_xlim([stepTime[1:].min(), stepTime.max()])

	print lineRna[:, 1].min(), lineRna[:, 1].max()
	
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
