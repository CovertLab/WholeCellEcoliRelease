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

NUMERICAL_ZERO = 1e-10 # used to cull very small numbers; could be chosen more rationally (i.e. minimum of one reaction per average cell size)

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	# Disabled
	return

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStep")
	reactionFluxes = fbaResults.readColumn("reactionFluxes")

	reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))

	fbaResults.close()

	# TODO: split figure output, perhaps using major clusterings

	REV_STR = " (reverse)"

	retainedReactionIndexes = []

	for reactionIndex, reactionID in enumerate(reactionIDs):
		if reactionID.endswith(REV_STR):
			forward_reactionID = reactionID.replace(REV_STR, "")

			forward_reactionIndex = np.where(forward_reactionID == reactionIDs)[0][0]

			reactionFluxes[:, forward_reactionIndex] -= reactionFluxes[:, reactionIndex]

		else:
			retainedReactionIndexes.append(reactionIndex)

	retainedReactionIndexes = np.array(retainedReactionIndexes)
	reactionFluxes = reactionFluxes[:, retainedReactionIndexes]
	reactionIDs = reactionIDs[retainedReactionIndexes]

	nonzero = (np.abs(reactionFluxes) >= NUMERICAL_ZERO).any(0)

	scaling = (np.mean(np.abs(reactionFluxes), 0) + 2 * np.std(np.abs(reactionFluxes), 0))

	# scaling = np.sqrt(np.corrcoef(reactionFluxes[:, nonzero].T)[np.identity(nonzero.sum(), np.bool)]) / reactionFluxes.shape[0]

	normalized = (
		reactionFluxes[:, nonzero]
		/ scaling[nonzero]
		).transpose()

	# with open(simDataFile, "rb") as f:
	# 	sim_data = cPickle.load(f)

	# idToName = {
	# 	reaction["id"]:reaction["name"]
	# 	for reaction in sim_data.process.metabolism.biochemicalReactions
	# 	}

	# reactionNames = np.array([
	# 	idToName[reactionID] for reactionID in reactionIDs
	# 	])

	# fig = plt.figure(figsize = (36, 48))
	# fig = plt.figure(figsize = (12, 16))

	linkage = sch.linkage(reactionFluxes[:, nonzero].T, metric = "correlation") #, method = "complete")
	linkage[:, 2] = np.fmax(linkage[:, 2], 0) # fixes rounding issues leading to negative distances

	sch.set_link_color_palette(['black'])

	fig = plt.figure(figsize = (72, 96))

	# grid = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0,width_ratios=[1,0.05])
	# dendro = sch.dendrogram(linkage, no_plot = True)

	grid = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0,width_ratios=[0.2, 1, 0.05])
	ax_dendro = fig.add_subplot(grid[0])

	dendro = sch.dendrogram(linkage, orientation="right", color_threshold = np.inf)

	ax_dendro.set_xticks([])
	ax_dendro.set_yticks([])
	ax_dendro.set_axis_off()

	index = dendro["leaves"]

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
	ax_mat.set_yticklabels(reactionIDs[nonzero][np.array(index)], size = 5)

	delta_t = time[1] - time[0]

	step_size = np.int64(5*60 / delta_t)

	xticks = np.arange(1, time.size, step_size)

	ax_mat.set_xticks(xticks)
	ax_mat.set_xticklabels(time[xticks-1]/60)

	ax_mat.set_xlabel("Time (min)")

	plt.title("Relative FBA reaction rates")

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

	plt.savefig(os.path.join(plotOutDir, plotOutFileName), dpi = 200)
	plt.close("all")


if __name__ == "__main__":
	# defaultSimDataFile = os.path.join(
	# 		wholecell.utils.constants.SERIALIZED_KB_DIR,
	# 		wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
	# 		)

	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	# parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])

	# Disabled until someone cares

	pass
