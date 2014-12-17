#!/usr/bin/env python
"""
Figure for the Simons Foundation statement of interest.  Combines other figures
and resizes for use in the manuscript.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/17/2014
"""

# TODO: separate the plotting code from the formatting code

import argparse
import os
import cPickle

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import gridspec

import wholecell.utils.constants

CATEGORICAL_COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[228,26,28],
	[55,126,184],
	[77,175,74],
	[152,78,163],
	[255,127,0],
	[255,255,51],
	[166,86,40],
	[247,129,191]
	]

CATEGORICAL_COLORS = [
	[colorValue/255. for colorValue in color]
	for color in CATEGORICAL_COLORS_256
	]

CMAP_COLORS_256 = [
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

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_256]

FIGURE_DIMENSIONS = (1.5*3, 1.5)

TICK_PAD = 2
LABEL_PAD = 2

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	plt.rc('font', size=3)
	plt.rc('xtick.major', pad=TICK_PAD)
	plt.rc('ytick.major', pad=TICK_PAD)

	fig = plt.figure(figsize = FIGURE_DIMENSIONS)

	grid = gridspec.GridSpec(2, 3)

	plotMassFractions(grid[:, 0], simOutDir, kbFile)

	plotRnaDistribution(grid[:, 1], simOutDir, kbFile)

	plotRnaAndProtein((grid[0, 2], grid[1, 2]), simOutDir, kbFile)

	#plt.show()

	grid.tight_layout(fig)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


def plotMassFractions(grids, simOutDir, kbFile):
	hdfFile = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))

	table = hdfFile.root.Mass

	# cell = table.col("cellMass")
	# cellDry = table.col("dryMass")
	protein = table.col("proteinMass")
	# rna = table.col("rnaMass")
	tRna = table.col("tRnaMass")
	rRna = table.col("rRnaMass")
	mRna = table.col("mRnaMass")
	dna = table.col("dnaMass")
	t = table.col("time")

	hdfFile.close()

	masses = np.vstack([
		protein/protein[0],
		rRna/rRna[0],
		tRna/tRna[0],
		mRna/mRna[0],
		dna/dna[0]
		]).T

	massLabels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA"]

	ax = plt.subplot(grids)

	ax.set_color_cycle(
		CATEGORICAL_COLORS[:masses.shape[1]]
		)

	plt.plot(t / 60., masses, linewidth = 2, alpha = 0.5)
	plt.plot(t / 60., masses, linewidth = 0.75, alpha = 1)
	plt.xlabel("Time (min)", labelpad = LABEL_PAD)
	plt.ylabel("Mass (normalized by t = 0 min)", labelpad = LABEL_PAD)
	plt.title("Biomass components")
	plt.axis([0, 60, 0.5, 2.5])

	plt.legend(massLabels, loc = "best")


def plotRnaDistribution(grids, simOutDir, kbFile):
	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	isMRna = kb.rnaData["isMRna"]

	rnaIds = kb.rnaData["id"][isMRna]

	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:

		names = bulkMoleculesFile.root.names
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules

		moleculeIds = names.moleculeIDs.read()

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]
	
	expectedCountsArbitrary = kb.rnaExpression['expression'][isMRna]

	expectedFrequency = expectedCountsArbitrary/expectedCountsArbitrary.sum()

	actualFrequency = 1.* rnaCountsBulk.T / rnaCountsBulk.sum(1)

	sorting = np.argsort(expectedFrequency)

	expectedFrequency = expectedFrequency[sorting]
	actualFrequency = actualFrequency[sorting, :]

	# Downsample for rendering purposes

	N_RNAS = 50
	N_TIMES = 20

	retainedRnas = np.int64(expectedFrequency.size - np.logspace(1, np.log10(expectedFrequency.size), N_RNAS))
	retainedTimes = np.arange(0, actualFrequency.shape[1], actualFrequency.shape[1] // N_TIMES)

	expectedFrequency = expectedFrequency[retainedRnas]
	actualFrequency = actualFrequency[retainedRnas, :][:, retainedTimes]

	maxValue = max(expectedFrequency.max(), actualFrequency.max()) * 1.1

	ax = plt.subplot(grids)

	plt.plot([0, maxValue], [0, maxValue], 'k--', linewidth = 0.2)

	# Plot style inspired by the prettyplotlib package (see:
	# http://blog.olgabotvinnik.com/post/58941062205/prettyplotlib-painlessly-create-beautiful-matplotlib)

	colors = np.array(
		(CATEGORICAL_COLORS*(expectedFrequency.size // len(CATEGORICAL_COLORS) + 1))[:expectedFrequency.size]
		).repeat(actualFrequency.shape[1], 0)

	plt.scatter(
		expectedFrequency.repeat(actualFrequency.shape[1]),
		actualFrequency.reshape(-1),
		c = colors,
		edgecolors = "none",
		#rasterized = True,
		alpha = 1,
		s = 1.5
		)

	plt.scatter(
		expectedFrequency, actualFrequency[:, -1],
		c = "none",
		edgecolors = np.array(
		(CATEGORICAL_COLORS*(expectedFrequency.size // len(CATEGORICAL_COLORS) + 1))[:expectedFrequency.size]
		),
		#linewidths = 0.5,
		#rasterized = True,
		#alpha = 0.5,
		marker = "+",
		#s = 2
		)

	plt.xlabel("Expected mRNA frequency", labelpad = LABEL_PAD)
	plt.ylabel("Actual mRNA frequency (multiple time steps)", labelpad = LABEL_PAD)
	plt.title("Simulation gene expression")

	plt.axis("scaled")
	plt.axis([0, maxValue, 0, maxValue])


def plotRnaAndProtein(grids, simOutDir, kbFile):
	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	rnaIds = kb.rnaData["id"][kb.rnaIndexToMonomerMapping]

	proteinIds = kb.monomerData["id"]

	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:

		names = bulkMoleculesFile.root.names
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules

		moleculeIds = names.moleculeIDs.read()

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]

		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

		proteinCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, proteinIndexes]

		time = bulkMolecules.read(0, None, 1, "time")

	time /= 60.

	# TODO: plots for all molecules (needs to be fast)
	# TODO: account for complex composition, other modifications

	index = 611
	name = "DNA polymerase"

	# TODO: functionalize this as a filled step plot

	stepTime = np.roll(np.repeat(time, 2), +1)[1:]

	stepRna = np.repeat(rnaCountsBulk[:, index], 2)[1:]
	stepProtein = np.repeat(proteinCountsBulk[:, index], 2)[1:]

	ax = plt.subplot(grids[0])

	polygonRna = np.vstack([
		np.hstack([stepTime[0], stepTime, stepTime[-1]]),
		np.hstack([0, stepRna, 0])
		]).T

	ax.add_patch(plt.Polygon(
		polygonRna,
		facecolor = CMAP_COLORS[3],
		edgecolor = "none"
		))

	lineRna = np.vstack([stepTime, stepRna]).T

	ax.add_patch(plt.Polygon(
		lineRna,
		facecolor = "none",
		edgecolor = CMAP_COLORS[1],
		closed = False,
		linewidth = 0.5
		))

	plt.axis("auto")

	#plt.xlabel("Time (min)")
	plt.ylabel("mRNA copies (number)", labelpad = LABEL_PAD)

	plt.title(name)

	ax = plt.subplot(grids[1])

	polygonProtein = np.vstack([
		np.hstack([stepTime[0], stepTime, stepTime[-1]]),
		np.hstack([0, stepProtein, 0])
		]).T

	ax.add_patch(plt.Polygon(
		polygonProtein,
		facecolor = CMAP_COLORS[-4],
		edgecolor = "none"
		))

	lineProtein = np.vstack([stepTime, stepProtein]).T

	ax.add_patch(plt.Polygon(
		lineProtein,
		facecolor = "none",
		edgecolor = CMAP_COLORS[-2],
		closed = False,
		linewidth = 0.5
		))

	plt.axis("auto")

	plt.xlabel("Time (min)", labelpad = LABEL_PAD)
	plt.ylabel("Protein copies (number)", labelpad = LABEL_PAD)


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
