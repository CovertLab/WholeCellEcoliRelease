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

HIGHLIGHT_PRE_POST_DNA_DIFF = False

fitterPredictionColor = "red"

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	return
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

	replicationResults = TableReader(os.path.join(simOutDir,"ReplicationData"))
	sequenceIdx = replicationResults.readColumn("sequenceIdx")
	replicationResults.close()

	activePolymerases = (sequenceIdx > -1).sum(axis=1)
	if np.cumsum(activePolymerases).any():
		replicationStartIdx = np.where(np.cumsum(activePolymerases) > 0)[0][0]
		replicationStopIdx = np.where(np.cumsum(activePolymerases) == np.cumsum(activePolymerases).max())[0][0]
	else:
		replicationStartIdx = 0
		replicationStopIdx = -1

	pointsToPlot, names = [], []
	for fluxName, predictedFlux in fitterPredictedFluxesDict.iteritems():
		reactionIdx = list(reactionIDs).index(fluxName)
		samplePoints = reactionFluxes[BURN_IN_PERIOD:, reactionIdx]
		if samplePoints.any():
			preReplicationPoints = samplePoints[:replicationStartIdx]
			midReplicationPoints = samplePoints[replicationStartIdx:replicationStopIdx]
			postReplicationPoints = samplePoints[replicationStopIdx:]
			samplePoints = [preReplicationPoints, midReplicationPoints, postReplicationPoints]
			pointsToPlot.append(samplePoints)
			names.append(fluxName)

	# Sort by mean flux
	pointsToPlot = np.array(sorted(pointsToPlot, key=lambda x: np.mean(np.abs(reduce(lambda z,y:list(z)+list(y), x)))))
	num_points = len(pointsToPlot)
	x_len = int(np.ceil(np.sqrt(num_points+4)))
	y_len = int(np.ceil((num_points+1)/x_len)+4)

	num_bins = 25

	fig = plt.figure(figsize=(np.ceil(2.5*x_len),np.ceil(2.5*y_len)))

	plt.suptitle("All Nonzero Reaction Fluxes After {} Step Burn-in, Sorted by Average Absolute Flux ({}) from Low to High.".format(BURN_IN_PERIOD, FLUX_UNITS.strUnit()), fontsize="xx-large")

	for idx, sample in enumerate(pointsToPlot):
		fluxName = names[idx]
		ax = plt.subplot(x_len,y_len,idx+1)
		plt.title(names[idx][:MAX_STRLEN],fontsize='xx-small')
		n, bins, patches = plt.hist(
			[x for x in sample],
			num_bins,
			histtype='bar',
			stacked=True,
			color=['green','blue', 'red'],
			label=['Before DNA replication', 'During DNA replication', ' After DNA replication'])
		n = reduce(lambda x,y:x+y, n)
		total_sample = np.array(reduce(lambda x,y:list(x)+list(y), sample))

		# Mann-Whitney U test if pre/during/post DNA replication come from the same distribution
		if HIGHLIGHT_PRE_POST_DNA_DIFF:
			u_value, p_value = stats.mannwhitneyu(sample[0], sample[1])
			if p_value < .05:
				for x in ax.spines:
					ax.spines[x].set_color('red')

		# Plot fitter-predicted flux on the axes
		ax.axvline(x=fitterPredictedFluxesDict[fluxName], color=fitterPredictionColor)

		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.tick_params(which = 'both', direction = 'out', labelsize=6)
		ax.set_xticks([])
		ymax = int(n.max())
		xmax = total_sample.max()
		xmin = total_sample.min()
		ax.set_ylim([0,ymax])
		ax.set_yticks([ymax])
		ax.set_yticklabels([str(ymax)])
		ax.set_xlim([xmin,xmax])
		ax.set_xticks([xmin,xmax])
		ax.set_xticklabels(["%.0e" % xmin,"%.0e" % xmax])


	ax = plt.subplot(x_len,y_len,idx+2)
	ax.xaxis.set_visible(False)
	ax.yaxis.set_visible(False)
	n, bins, patches = plt.hist(
			[[0] for x in sample],
			num_bins,
			histtype='bar',
			stacked=True,
			color=['green','blue', 'red'],
			label=['Before DNA replication', 'During DNA replication', ' After DNA replication'])
	plt.legend(loc='center left', fontsize="xx-large")

	ax = plt.subplot(x_len,y_len,idx+5)
	ax.axis("off")
	# ax.xaxis.set_visible(False)
	# ax.yaxis.set_visible(False)
	plt.text(0,.8, "{} vertical line is fitter predicted value.".format(fitterPredictionColor), fontsize="xx-large", color=fitterPredictionColor)

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