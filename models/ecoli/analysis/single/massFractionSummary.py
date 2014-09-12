#!/usr/bin/env python
"""
Plot mass fractions

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

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

	h = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))

	table = h.root.Mass

	# cell = table.col("cellMass")
	# cellDry = table.col("dryMass")
	protein = table.col("proteinMass")
	# rna = table.col("rnaMass")
	tRna = table.col("tRnaMass")
	rRna = table.col("rRnaMass")
	mRna = table.col("mRnaMass")
	dna = table.col("dnaMass")
	t = table.col("time")

	masses = np.vstack([
		protein/protein[0],
		rRna/rRna[0],
		tRna/tRna[0],
		mRna/mRna[0],
		dna/dna[0]
		]).T

	massLabels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA"]

	plt.figure(figsize = (8.5, 11))

	plt.rc('axes', color_cycle=COLORS)

	plt.plot(t / 60., masses, linewidth = 2)
	plt.xlabel("Time (min)")
	plt.ylabel("Mass (normalized by t = 0 min)")
	plt.title("Biomass components")
	plt.axis([0, 60, 0.5, 2.5])

	plt.legend(massLabels, loc = "best")

	# plt.show()

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

	h.close()

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
