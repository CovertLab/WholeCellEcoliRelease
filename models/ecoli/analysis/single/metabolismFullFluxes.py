#!/usr/bin/env python
"""
@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/04/2016
"""

from __future__ import division

import argparse
import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
import scipy.cluster.hierarchy as sch

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

CMAP_COLORS_255 = [
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	reactionFluxes = fbaResults.readColumn("reactionFluxes")

	reactionIDs = np.array(fbaResults.readAttribute('reactionIDs'))

	fbaResults.close()

	fig = plt.figure(figsize = (80, 40))

	grid = gridspec.GridSpec(1,2,wspace=0.0,hspace=0.0,width_ratios=[1,0.1])

	normalized = (
		reactionFluxes
		/ (np.mean(np.abs(reactionFluxes), 0) + 2 * np.std(np.abs(reactionFluxes), 0))
		).transpose()

	linkage = sch.linkage(reactionFluxes.T)
	linkage[:, 2] = np.fmax(linkage[:, 2], 0) # fixes rounding issues leading to negative distances

	sch.set_link_color_palette(['black'])

	index = sch.leaves_list(linkage)

	ax_mat = fig.add_subplot(grid[0])

	cmap = colors.LinearSegmentedColormap.from_list(
		"white to blue with upper extreme",
		CMAP_COLORS
		)

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
	ax_mat.set_yticklabels(reactionIDs[index], size = 1)

	xticks = np.linspace(0,time.size-1, 5, dtype=np.int)
	ax_mat.set_xticks(xticks)
	ax_mat.set_xticklabels(np.round(time[xticks]/60.), size=40)
	ax_mat.set_xlabel("Time (min)", size=40)

	plt.title("Full network FBA reaction fluxes", size =40)

	ax_cmap = fig.add_subplot(grid[1])

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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
