#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

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

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	mass = TableReader(os.path.join(simOutDir, "Mass"))

	# cell = mass.readColumn("cellMass")
	# cellDry = mass.readColumn("dryMass")
	protein = mass.readColumn("proteinMass")
	# rna = mass.readColumn("rnaMass")
	tRna = mass.readColumn("tRnaMass")
	rRna = mass.readColumn("rRnaMass")
	mRna = mass.readColumn("mRnaMass")
	dna = mass.readColumn("dnaMass")

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	f = plt.figure(figsize = (1.25, 0.78), frameon = False)
	ax = f.add_axes([0, 0, 1, 1])
	ax.axis("off")

	ax.set_color_cycle(COLORS)

	ax.plot(t / 60., mRna, linewidth = 2)
	ax.set_ylim([mRna.min(), mRna.max()])

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
