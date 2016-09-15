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
from wholecell.utils import units

import mpld3
from mpld3 import plugins, utils

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE, plotSplom

NUMERICAL_ZERO = 1e-15

BURN_IN_PERIOD = 150

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	sim_data = cPickle.load(open(simDataFile, "rb"))
	
	massListener = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = massListener.readColumn("cellMass") * units.fg
	dryMass = massListener.readColumn("dryMass") * units.fg
	massListener.close()

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
	reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
	fbaResults.close()
	
	cellDensity = sim_data.constants.cellDensity
	dryMassFracAverage = np.mean(dryMass / cellMass)

	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])
	toya_fluxes_dict = dict(zip(toya_reactions, toya_fluxes))

	# Clip all values less than numerical zero to zero
	reactionFluxes[reactionFluxes < NUMERICAL_ZERO] = 0

	fitterPredictedFluxesDict = {key:value.asNumber(FLUX_UNITS) for key, value in sim_data.process.metabolism.predictedFluxesDict.iteritems() if key in reactionIDs}
	for key, value in fitterPredictedFluxesDict.iteritems():
		if np.abs(value) < NUMERICAL_ZERO:
			fitterPredictedFluxesDict[key] = 0

	scatterArrayActual, scatterArrayPredicted, scatterArrayActualStd, toyaObservedFluxes, labels = [], [], [], [], []
	for fluxName, toyaFlux in toya_fluxes_dict.iteritems():
		reactionIdx = list(reactionIDs).index(fluxName)
		samplePoints = reactionFluxes[BURN_IN_PERIOD:, reactionIdx]
		scatterArrayActual.append(np.mean(samplePoints))
		scatterArrayPredicted.append(fitterPredictedFluxesDict[fluxName])
		scatterArrayActualStd.append(np.std(samplePoints))
		toyaObservedFluxes.append(toyaFlux.asNumber(FLUX_UNITS))
		labels.append(fluxName)
	
	scatterArrayActual = np.array(scatterArrayActual)
	scatterArrayPredicted = np.array(scatterArrayPredicted)
	toyaObservedFluxes = np.array(toyaObservedFluxes)

	arrayOfdataArrays = [scatterArrayActual, scatterArrayPredicted, toyaObservedFluxes]
	arrayOfdataStdArrays = [scatterArrayActualStd, None, None]

	names = ["WCM Flux {}".format(FLUX_UNITS.strUnit()), "Fitter Prediction {}".format(FLUX_UNITS.strUnit()), "Toya et al Measurement {}".format(FLUX_UNITS.strUnit())]

	fig = plt.figure(figsize=(30,30))
	plt.suptitle("Actual vs. Fitter Predicted vs. Toya Observed Fluxes, {} Step Burn-in.".format(BURN_IN_PERIOD))

	fig = plotSplom(arrayOfdataArrays, nameArray=names, stdArrays=arrayOfdataStdArrays, labels=labels, fig=fig, plotCorrCoef=True, htmlPlot=True)

	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)

	# Error bars don't work for the HTML plot, so save it with just points, then add error bars.
	fig = plotSplom(arrayOfdataArrays, nameArray=names, stdArrays=arrayOfdataStdArrays, fig=fig, plotCorrCoef=True)

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
