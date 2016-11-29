#!/usr/bin/env python
"""
Plot the components of the objective associated with homeostatic targets.

@date: Created 9/28/2016
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

from wholecell.utils import units
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

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

from wholecell.analysis.plotting_tools import CMAP_COLORS_255
CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_UNDER = [1, 0.2, 0.75]
CMAP_OVER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaData = TableReader(os.path.join(simOutDir, "FBAResults"))
	
	kineticTargetFluxNames = np.array(fbaData.readAttribute("kineticTargetFluxNames"))
	kineticObjectiveComponents = fbaData.readColumn("kineticObjectiveValues")
	homeostaticObjectiveWeights = fbaData.readColumn("homeostaticObjectiveWeight")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fbaData.close()

	fig = plt.figure(figsize = (12, 6))

	grid = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0,width_ratios=[0.25,1,0.05])

	ax_dendro = fig.add_subplot(grid[0])

	normalized = (
		kineticObjectiveComponents
		/ (np.mean(np.abs(kineticObjectiveComponents), 0) + 2 * np.std(np.abs(kineticObjectiveComponents), 0))
		).transpose()

	linkage = sch.linkage(kineticObjectiveComponents.T, metric = "correlation")
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
	ax_mat.set_yticklabels(kineticTargetFluxNames[np.array(index)], size = 4)

	xticks = np.arange(time[0], time[-1], 720)
	tickInd = np.zeros(len(xticks))
	for i, tick in enumerate(xticks):
		for j, t in enumerate(time):
			if t >= tick:
				tickInd[i] = j
				break

	xtickmin = np.arange(0, time[-1] / 3600., 0.2)

	ax_mat.set_xticks(tickInd)
	ax_mat.set_xticklabels(xtickmin, size = 4)

	ax_mat.set_xlabel("Time (hr)", size = 6)

	plt.title("Kinetic objective component values", size = 8)

	ax_cmap = fig.add_subplot(grid[2])

	gradient = np.array((np.arange(0, 100)/100).tolist() + [+2,]*5, ndmin=2).transpose()

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