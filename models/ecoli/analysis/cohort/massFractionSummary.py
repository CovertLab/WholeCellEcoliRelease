#!/usr/bin/env python

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

	simOutSubDirs = sorted([
		os.path.join(simOutDir, item, "simOut")
		for item in os.listdir(simOutDir)
		if os.path.isdir(os.path.join(simOutDir, item)) and item not in {"kb", "metadata"}
		])

	time = None

	for simIndex, simOutSubDir in enumerate(simOutSubDirs):

		with tables.open_file(os.path.join(simOutSubDir, "Mass.hdf")) as hdfFile:

			table = hdfFile.root.Mass

			if time is None:
				time = table.col("time")

				shape = (time.size, len(simOutSubDirs))

				protein = np.empty(shape, order = "F")
				tRna = np.empty(shape, order = "F")
				rRna = np.empty(shape, order = "F")
				mRna = np.empty(shape, order = "F")
				dna = np.empty(shape, order = "F")

			mRna[:, simIndex] = table.col("mRnaMass")
			tRna[:, simIndex] = table.col("tRnaMass")
			rRna[:, simIndex] = table.col("rRnaMass")
			protein[:, simIndex] = table.col("proteinMass")
			dna[:, simIndex] = table.col("dnaMass")

	massesMean = np.vstack([
		mRna.mean(1),
		tRna.mean(1),
		rRna.mean(1),
		protein.mean(1),
		dna.mean(1)
		]).T

	# massesStd = np.vstack([
	# 	mRna.std(1),
	# 	tRna.std(1),
	# 	rRna.std(1),
	# 	protein.std(1),
	# 	dna.std(1)
	# 	]).T

	massesMax = np.vstack([
		mRna.max(1),
		tRna.max(1),
		rRna.max(1),
		protein.max(1),
		dna.max(1)
		]).T

	massesMin = np.vstack([
		mRna.min(1),
		tRna.min(1),
		rRna.min(1),
		protein.min(1),
		dna.min(1)
		]).T

	time /= 60

	# Normalize

	normalization = massesMean[0, :].copy()

	massesMean /= normalization
	massesMin /= normalization
	massesMax /= normalization
	# massesStd /= normalization

	massLabels = ["mRNA", "tRNA", "rRNA", "Protein", "DNA"]

	plt.figure(figsize = (8.5, 11))

	# plt.rc('axes', color_cycle=COLORS)

	for i in xrange(len(massLabels)):
		plt.fill_between(time, massesMin[:, i], massesMax[:, i],
			alpha = 0.25, facecolor = COLORS[i], edgecolor = "none")

		# plt.fill_between(time, massesMean[:, i] - massesStd[:, i], massesMean[:, i] + massesStd[:, i],
		# 	alpha = 0.5, facecolor = COLORS[i], edgecolor = "none")

	plt.gca().set_color_cycle(COLORS[:len(massLabels)])

	plt.plot(time, massesMean, linewidth = 2)
	plt.xlabel("Time (min)")
	plt.ylabel("Mass (normalized by t = 0 min)")
	plt.title("Biomass components")
	plt.axis([0, 60, 0.5, 2.5])

	plt.legend(massLabels, loc = "best")

	# plt.show()

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

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
