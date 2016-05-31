#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2016
"""
from __future__ import division

import argparse
import os
import sys
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance

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

CMAP_COLORS = [[shade/255. for shade in color] for color in reversed(CMAP_COLORS_255)]
CMAP_OVER = [1, 0.2, 0.75]
CMAP_UNDER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# sys.setrecursionlimit(3500)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	dualValues = np.array(fbaResults.readColumn("dualValues"))
	moleculeIDs = np.array(fbaResults.readAttribute("outputMoleculeIDs"))
	fbaResults.close()

	fig = plt.figure(figsize = (30, 15))

	grid = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0,width_ratios=[0.25,1,0.1])

	ax_dendro = fig.add_subplot(grid[0])

	scaling = (np.mean(np.abs(dualValues), 0) + 2 * np.std(np.abs(dualValues), 0))

	nonzero = (scaling != 0)

	normalized = (
		dualValues[:, nonzero]
		/ scaling[nonzero]
		).transpose()

	linkage = sch.linkage(dualValues[:, nonzero].T, metric = "correlation")
	linkage[:, 2] = np.fmax(linkage[:, 2], 0) # fixes rounding issues leading to negative distances

	sch.set_link_color_palette(['black'])

	dendro = sch.dendrogram(linkage, orientation="right", color_threshold = np.inf)
	index = dendro["leaves"]

	ax_dendro.set_xticks([])
	ax_dendro.set_yticks([])
	ax_dendro.set_axis_off()

	ax_mat = fig.add_subplot(grid[1])

	cmap = colors.LinearSegmentedColormap.from_list(
		"red to blue with extremes",
		CMAP_COLORS
		)

	cmap.set_under(CMAP_UNDER)
	cmap.set_over(CMAP_OVER)

	norm = colors.Normalize(vmin = -1, vmax = +1)

	ax_mat.imshow(
		normalized[index, :],
		aspect = "auto",
		interpolation='nearest',
		origin = "lower",
		cmap = cmap,
		norm = norm
		)

	ax_mat.set_yticks(np.arange(len(index)))
	ax_mat.set_yticklabels(moleculeIDs[np.array(index)], size = 5)

	xticks = np.linspace(0,time.size-1, 5, dtype=np.int)
	ax_mat.set_xticks(xticks)
	ax_mat.set_xticklabels(np.round(time[xticks]/60.,decimals=1))
	ax_mat.set_xlabel("Time (min)")

	plt.title("FBA Reduced Prices (Red = wants less, Blue = wants more)")

	ax_cmap = fig.add_subplot(grid[2])

	gradient = np.array([-2,]*5 + (np.arange(-100, 100)/100).tolist() + [+2,]*5, ndmin=2).transpose()

	ax_cmap.imshow(
		gradient,
		aspect = "auto",
		interpolation = "nearest",
		origin = "lower",
		cmap = cmap,
		norm = norm
		)

	ax_cmap.set_xticks([])
	ax_cmap.set_yticks([])

	grid.tight_layout(fig)

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
	