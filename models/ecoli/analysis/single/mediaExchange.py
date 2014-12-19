#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import argparse
import os
import cPickle

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance

import wholecell.utils.constants

FLUX_UNITS = "M/s"

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
CMAP_UNDER = [1, 0.2, 0.75]
CMAP_OVER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	with tables.open_file(os.path.join(simOutDir, "FBAResults.hdf")) as h5file:
		time = h5file.root.FBAResults.col("time")
		timeStep = h5file.root.FBAResults.col("timeStep")
		externalExchangeFluxes = h5file.root.FBAResults.col("externalExchangeFluxes")

		names = h5file.root.names
		externalMoleculeIDs = np.array(names.externalMoleculeIDs.read())

	fig = plt.figure(figsize = (30, 15))

	grid = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0,width_ratios=[0.25,1,0.1])

	ax_dendro = fig.add_subplot(grid[0])

	scaling = (np.mean(np.abs(externalExchangeFluxes), 0) + 2 * np.std(np.abs(externalExchangeFluxes), 0))

	nonzero = (scaling != 0)

	normalized = (
		externalExchangeFluxes[:, nonzero]
		/ scaling[nonzero]
		).transpose()

	linkage = sch.linkage(externalExchangeFluxes[:, nonzero].T, metric = "correlation")
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
	ax_mat.set_yticklabels(externalMoleculeIDs[nonzero][np.array(index)], size = 5)

	delta_t = time[1] - time[0]

	step_size = np.int64(5*60 / delta_t)

	xticks = np.arange(1, time.size, step_size)

	ax_mat.set_xticks(xticks)
	ax_mat.set_xticklabels(time[xticks-1]/60)

	ax_mat.set_xlabel("Time (min)")

	plt.title("Relative FBA media utilization rates (red = secretion, blue = uptake)")

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
