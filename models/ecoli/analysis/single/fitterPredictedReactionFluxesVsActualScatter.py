#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/14/2016
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

	# Clip all values less than numerica zero to zero
	reactionFluxes[reactionFluxes < NUMERICAL_ZERO] = 0

	fitterPredictedFluxesDict = {key:value.asNumber(FLUX_UNITS) for key, value in sim_data.process.metabolism.predictedFluxesDict.iteritems() if key in reactionIDs}
	for key, value in fitterPredictedFluxesDict.iteritems():
		if np.abs(value) < NUMERICAL_ZERO:
			fitterPredictedFluxesDict[key] = 0

	scatterArrayActual, scatterArrayPredicted, scatterArrayPredictedStd, labels = [], [], [], []
	for fluxName, predictedFlux in fitterPredictedFluxesDict.iteritems():
		reactionIdx = list(reactionIDs).index(fluxName)
		samplePoints = reactionFluxes[BURN_IN_PERIOD:, reactionIdx]
		scatterArrayActual.append(predictedFlux)
		scatterArrayPredicted.append(np.mean(samplePoints))
		scatterArrayPredictedStd.append(np.std(samplePoints))
		labels.append(fluxName)
	
	scatterArrayActual = np.array(scatterArrayActual)
	scatterArrayPredicted = np.array(scatterArrayPredicted)
	scatterArrayPredictedStd = np.array(scatterArrayPredictedStd)
	correlationCoefficient = np.corrcoef(scatterArrayActual, scatterArrayPredicted)[0,1]
	logCorrelationCoefficient = np.corrcoef(np.log10(scatterArrayActual + 1), np.log10(scatterArrayPredicted + 1))[0,1]

	fig = plt.figure(figsize=(7.5,11))

	plt.suptitle("Fitter Predicted Fluxes, {} Step Burn-in.".format(BURN_IN_PERIOD))

	plt.subplot(2,1,1)
	plt.title("Pearson R = {:.2}".format(correlationCoefficient))
	points = plt.scatter(scatterArrayPredicted, scatterArrayActual)
	plt.xlabel("Predicted Flux {}".format(FLUX_UNITS.strUnit()))
	plt.ylabel("Actual Flux {}".format(FLUX_UNITS.strUnit()))

	tooltip = plugins.PointLabelTooltip(points, labels)
	plugins.connect(fig, tooltip)

	plt.subplot(2,1,2)
	plt.title("Pearson R = {:.2}".format(logCorrelationCoefficient))
	points = plt.scatter(np.log10(scatterArrayPredicted), np.log10(scatterArrayActual))
	plt.xlabel("Log10 Predicted Flux {}".format(FLUX_UNITS.strUnit()))
	plt.ylabel("Log10 Actual Flux {}".format(FLUX_UNITS.strUnit()))

	tooltip = plugins.PointLabelTooltip(points, labels)
	plugins.connect(fig, tooltip)

	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)

	# Error bars don't work for the HTML plot, so save it with just points, then add error bars.
	plt.subplot(2,1,1)
	points = plt.errorbar(scatterArrayPredicted, scatterArrayActual, yerr=scatterArrayPredictedStd, fmt='o')

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
