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

NUMERICAL_ZERO = 1e-15

BURN_IN_PERIOD = 150

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

	p_values, nonzero_p_values = [], []
	for fluxName, predictedFlux in fitterPredictedFluxesDict.iteritems():
		reactionIdx = list(reactionIDs).index(fluxName)
		samplePoints = reactionFluxes[BURN_IN_PERIOD:, reactionIdx]
		degsOfFreedom = len(samplePoints) - 1
		if samplePoints.any():
			t_value = (np.mean(samplePoints) - predictedFlux) / (np.std(samplePoints) / np.sqrt(degsOfFreedom + 1))
			p_value = stats.t.sf(np.abs(t_value), degsOfFreedom)*2.
			p_values.append(p_value)
			nonzero_p_values.append(p_value)
		else:
			t_value = (np.mean(samplePoints) - predictedFlux) / (NUMERICAL_ZERO / np.sqrt(degsOfFreedom + 1))
			p_value = stats.t.sf(np.abs(t_value), degsOfFreedom)*2.
			p_values.append(p_value)

	p_values = np.array(p_values)
	nonzero_p_values = np.array(nonzero_p_values)

	fig = plt.figure(figsize=(7.5,11))

	plt.suptitle("T-test of Predicted Flux vs Observed Fluxes, {} Step Burn-in.".format(BURN_IN_PERIOD))
	
	ax1 = plt.subplot(2,1,1)
	plt.title("All Fluxes ({})".format(len(p_values)))
	n, bins, patches = plt.hist(p_values, 20, range=(0,1), normed=False)
	plt.ylabel("Number of Fluxes")
	percentOver05 = (p_values > .05).mean()
	plt.text(.5, .95*(n.max()), "{0:.1f}% over .05".format(percentOver05*100.), horizontalalignment='center')

	ax2 = plt.subplot(2,1,2, sharex=ax1)
	plt.title("Nonzero Fluxes Only ({} WCM fluxes always zero)".format(len(p_values) - len(nonzero_p_values)))
	n, bins, patches = plt.hist(nonzero_p_values, 20, range=(0,1), normed=False)
	plt.xlabel("P-value")
	plt.ylabel("Number of Fluxes")
	percentOver05 = (nonzero_p_values > .05).mean()
	plt.text(.5, .95*(n.max()), "{0:.1f}% over .05".format(percentOver05*100.), horizontalalignment='center')

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
