#!/usr/bin/env python
"""
Cluster and plot counts for enzymes which catalyze metabolic reactions

@date: Created 5/18/2016
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

VERBOSE = True

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

	enzymeCountsInitRaw = enzymeKineticsdata.readColumn("enzymeCountsInit")
	countsToMolarArray = enzymeKineticsdata.readColumn("countsToMolar")
	reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")

	# Read time info from the listener
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))

	enzymeNames = np.array(sim_data.process.metabolism.enzymeNames)
	reactionEnzymes = sim_data.process.metabolism.reactionEnzymes
	catalyticEnzymes = set([enzymeName for enzymeList in reactionEnzymes.values() for enzymeName in enzymeList])
	
	catalyticEnzymeNames = []
	indices = []
	for enzymeName in enzymeNames:
		if enzymeName in catalyticEnzymes:
			catalyticEnzymeNames.append(enzymeName)
			index = list(enzymeNames).index(enzymeName)
			indices.append(index)
	catalyticEnzymeNames = np.array(catalyticEnzymeNames)		

	catalyticEnzymeCounts = enzymeCountsInitRaw[:,indices]

	alwaysZero = enzymeNames[np.where(np.sum(catalyticEnzymeCounts, axis=0) == 0)]

	catalyticEnzymeConcentrations = catalyticEnzymeCounts * countsToMolarArray.reshape(-1,1)

	if VERBOSE: print "{} catalytic enzymes never produced: {}".format(len(alwaysZero), alwaysZero)


	fig = plt.figure(figsize = (120, 60))

	grid = gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0,width_ratios=[0.25,1,0.1])

	ax_dendro = fig.add_subplot(grid[0])

	normalized = (
		catalyticEnzymeConcentrations
		/ (np.mean(np.abs(catalyticEnzymeConcentrations), 0) + 2 * np.std(np.abs(catalyticEnzymeConcentrations), 0))
		).transpose()

	linkage = sch.linkage(catalyticEnzymeConcentrations.T, metric = "correlation")
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
	ax_mat.set_yticklabels(catalyticEnzymeNames[np.array(index)], size = 5)

	delta_t = time[1] - time[0]

	step_size = np.int64(5*60 / delta_t)

	xticks = np.arange(1, time.size, step_size)

	ax_mat.set_xticks(xticks)
	ax_mat.set_xticklabels(time[xticks-1]/60)

	ax_mat.set_xlabel("Time (min)")

	plt.title("Relative counts of catalytic enzymes")

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