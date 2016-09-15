#!/usr/bin/env python
"""
Plot distance from target value for all FBA fluxes

@date: Created 7/29/2016
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os
import itertools
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import FLUX_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE

NUMERICAL_ZERO = 1e-20

MAX_STR_LEN = 30
MIN_RELATIVE_ERROR_TO_SHOW = 1e-5
MOVING_AVE_WINDOW_SIZE = 100
BURN_IN_PERIOD = 150
NUM_ERRORS_TO_SHOW = 5

TARGET_RELATIVE_DIFFERENCE = 1

PERCENT_ON_TARGET_THRESHOLD = 1.0

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	errorRelativeDifferences = enzymeKineticsdata.readColumn("kineticTargetRelativeDifferences")
	errorFluxNames = enzymeKineticsdata.readAttribute("kineticTargetFluxNames")
	oneSidedErrorFluxNames = enzymeKineticsdata.readAttribute("kineticOneSidedTargets")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	enzymeKineticsdata.close()

	# Build a mapping from reaction to color
	idToColor = {}
	for errorFluxID, color in itertools.izip(errorFluxNames, itertools.cycle(COLORS_LARGE)):
		idToColor[errorFluxID] = color

	# Determine where to set the threshold to hit the target number of most deviant fluxes
	with warnings.catch_warnings() as w:
		warnings.simplefilter("ignore")
		# Smooth the error fluxes with a moving average
		smoothedErrorFluxes = np.zeros_like(errorRelativeDifferences[BURN_IN_PERIOD:])
		for idx, timeCourse in enumerate(errorRelativeDifferences[BURN_IN_PERIOD:].T):
			smoothedTimeCourse = np.convolve(timeCourse, np.ones((MOVING_AVE_WINDOW_SIZE,))/MOVING_AVE_WINDOW_SIZE, mode='same')
			smoothedErrorFluxes[:,idx] = smoothedTimeCourse
		maxima = np.nanmax(np.abs(smoothedErrorFluxes), axis=0)
	maxima = np.delete(maxima, np.where(np.isnan(maxima)))
	threshold_value = sorted(maxima)[-(NUM_ERRORS_TO_SHOW+1)]

	# Track fluxes with errors always below numerical zero
	zeroErrorFluxes = set()

	# Track fluxes within a specified range of the target
	withinTargetErrorFluxes = set()

	plt.figure(figsize = (17, 11))
	plt.suptitle("Relative Difference from Target Flux for Reactions with Kinetic Estimates")

	for idx, (errorFluxID, errorFluxTimeCourse) in enumerate(zip(errorFluxNames, errorRelativeDifferences.T)):

		if (np.abs(errorFluxTimeCourse[BURN_IN_PERIOD:]) <= NUMERICAL_ZERO).mean() >= PERCENT_ON_TARGET_THRESHOLD:
			zeroErrorFluxes.add(errorFluxID)
			withinTargetErrorFluxes.add(errorFluxID)
			continue

		if (np.abs(errorFluxTimeCourse[BURN_IN_PERIOD:]) <= TARGET_RELATIVE_DIFFERENCE).mean() >= PERCENT_ON_TARGET_THRESHOLD:
			withinTargetErrorFluxes.add(errorFluxID)

		plt.subplot(2,2,1)
		plt.plot(time[BURN_IN_PERIOD:] / 60., errorFluxTimeCourse[BURN_IN_PERIOD:], label=errorFluxNames[idx][:MAX_STR_LEN], color=idToColor[errorFluxID])
	
		plt.subplot(2,2,3)
		plt.plot(time[BURN_IN_PERIOD:] / 60., np.log10(np.abs(errorFluxTimeCourse[BURN_IN_PERIOD:])), label=errorFluxNames[idx][:MAX_STR_LEN], color=idToColor[errorFluxID])
		plt.xlabel("Time (min)")
		plt.ylabel("Log10 Absolute Flux Error ({})".format(FLUX_UNITS.strUnit()))

		if ((np.abs(smoothedErrorFluxes[:,idx]) >= threshold_value).any()) and np.abs(errorFluxTimeCourse[BURN_IN_PERIOD:]).max() > MIN_RELATIVE_ERROR_TO_SHOW:
			plt.subplot(2,2,2)
			plt.plot(time[BURN_IN_PERIOD:] / 60., errorFluxTimeCourse[BURN_IN_PERIOD:], label=errorFluxNames[idx][:MAX_STR_LEN], color=idToColor[errorFluxID])

			plt.subplot(2,2,4)
			plt.plot(time[BURN_IN_PERIOD:] / 60., np.log10(np.abs(errorFluxTimeCourse[BURN_IN_PERIOD:])), label=errorFluxNames[idx][:MAX_STR_LEN], color=idToColor[errorFluxID])
			plt.xlabel("Time (min)")

	plt.subplot(2,2,1)
	plt.ylabel("Relative Difference From Target")

	plt.subplot(2,2,2)
	plt.title("{} max absolute relative diff after moving mean smoothing window {}, burn-in {}. Must be >{}.".format(NUM_ERRORS_TO_SHOW, MOVING_AVE_WINDOW_SIZE, BURN_IN_PERIOD, MIN_RELATIVE_ERROR_TO_SHOW),fontsize='x-small')

	plt.subplot(2,2,3)
	plt.title("{} fluxes ({:.1f}%) with relative difference <={} {}% of the time".format(len(withinTargetErrorFluxes), 100*(len(withinTargetErrorFluxes)/len(errorFluxNames)), TARGET_RELATIVE_DIFFERENCE,100*PERCENT_ON_TARGET_THRESHOLD), fontsize='small')
	plt.xlabel("Time (min)")
	plt.ylabel("Log10 Absolute Relative Difference")

	plt.subplot(2,2,4)
	plt.xlabel("Time (min)")
	plt.legend(fontsize="xx-small", framealpha=.5, loc="best")

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
