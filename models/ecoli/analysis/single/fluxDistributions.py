#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/20/2016
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import itertools
import warnings
from scipy import stats

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

import mpld3
from mpld3 import plugins, utils

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE

BURN_IN_PERIOD = 150
MAX_STRLEN = 20

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))
	
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
	reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
	fbaResults.close()

	fitterPredictedFluxesDict = {key:value.asNumber(FLUX_UNITS) for key, value in sim_data.process.metabolism.predictedFluxesDict.iteritems() if key in reactionIDs}

	pointsToPlot, names = [], []
	for fluxName, predictedFlux in fitterPredictedFluxesDict.iteritems():
		reactionIdx = list(reactionIDs).index(fluxName)
		samplePoints = reactionFluxes[BURN_IN_PERIOD:, reactionIdx]
		if samplePoints.any():
			pointsToPlot.append(samplePoints)
			names.append(fluxName)

	fig = plt.figure(figsize=(50,50))

	num_points = len(pointsToPlot)
	x_len = int(np.ceil(np.sqrt(num_points)))
	y_len = int(np.ceil(num_points/x_len) + 1)

	num_bins=30

	plt.suptitle("All nonzero reaction fluxes", fontsize="xx-large")

	for idx in xrange(1,num_points-1):
		ax = plt.subplot(x_len,y_len,idx)
		plt.title(names[idx-1][:MAX_STRLEN],fontsize='xx-small')
		n, bins, patches = plt.hist(pointsToPlot[idx-1],num_bins)
		ax.xaxis.set_visible(False)
	
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=6)
		ax.set_xticks([])
		ymax = int(n.max())
		ax.set_ylim([0,ymax])
		ax.set_yticks([ymax])
		ax.set_yticklabels([str(ymax)])

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
