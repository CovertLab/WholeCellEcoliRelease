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

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE

BURN_IN_PERIOD = 150

NUMERICAL_ZERO = 1e-15

MIN_PERCENT_DEVIATION_TO_SHOW = .001
MAX_NUM_DEVIANTS_TO_SHOW = 5

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	idToColor = {}

	sim_data = cPickle.load(open(simDataFile, "rb"))
	
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
	reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
	fbaResults.close()

	fitterPredictedFluxesDict = sim_data.process.metabolism.predictedFluxesDict
	fitterPredictedFluxes = np.array([fitterPredictedFluxesDict[reactionID].asNumber(FLUX_UNITS) for reactionID in reactionIDs])

	# Clip reaction fluxes which are less than numerical zero to numerical zero
	reactionFluxes[np.abs(reactionFluxes) < NUMERICAL_ZERO] = 0
	fitterPredictedFluxes[np.abs(fitterPredictedFluxes) < NUMERICAL_ZERO] = 0

	percentDeviations = (fitterPredictedFluxes - reactionFluxes) / fitterPredictedFluxes

	# Don't plot differences for fluxes for which the prediction was zero
	percentDeviations[np.isinf(percentDeviations)] = np.nan
	
	# Build a mapping from reaction to color
	idToColor = {}
	for reactionID, color in itertools.izip(reactionIDs, itertools.cycle(COLORS_LARGE)):
		idToColor[reactionID] = color
	
	plt.figure(figsize = (17, 11))

	# Determine where to set the threshold to hit the target number of most deviant fluxes
	with warnings.catch_warnings() as w:
		warnings.simplefilter("ignore")
		maxima = np.nanmax(np.abs(percentDeviations[BURN_IN_PERIOD:]), axis=0)
	maxima = np.delete(maxima, np.where(np.isnan(maxima)))
	threshold_value = sorted(maxima)[-(MAX_NUM_DEVIANTS_TO_SHOW+1)]

	for idx, (reactionID, percentDeviation) in enumerate(zip(reactionIDs, percentDeviations.T)):
		# Unadjusted
		plt.subplot(2,2,1)
		plt.plot(time[BURN_IN_PERIOD:]/60., percentDeviation[BURN_IN_PERIOD:], label=reactionID, color=idToColor[reactionID])

		# Log scale
		plt.subplot(2,2,3)
		plt.plot(time[BURN_IN_PERIOD:]/60., np.log10(np.abs(percentDeviation[BURN_IN_PERIOD:])), label=reactionID, color=idToColor[reactionID])

		if ((np.abs(percentDeviation[BURN_IN_PERIOD:]) > threshold_value).any()) and np.abs(percentDeviation[BURN_IN_PERIOD:]).max() > MIN_PERCENT_DEVIATION_TO_SHOW:
			# Unadjusted
			plt.subplot(2,2,2)
			plt.plot(time[BURN_IN_PERIOD:]/60., percentDeviation[BURN_IN_PERIOD:], label=reactionID, color=idToColor[reactionID])

			# Log scale
			plt.subplot(2,2,4)
			plt.plot(time[BURN_IN_PERIOD:]/60., np.log10(np.abs(percentDeviation[BURN_IN_PERIOD:])), label=reactionID, color=idToColor[reactionID])

	plt.suptitle("Reaction Fluxes")
	plt.subplot(2,2,1)
	plt.ylabel('Relative Difference from Prediction')
	plt.subplot(2,2,2)
	plt.title("Top {} fluxes by max percent deviation from prediction after {} step burn-in. Must be over {}.".format(MAX_NUM_DEVIANTS_TO_SHOW, BURN_IN_PERIOD, MIN_PERCENT_DEVIATION_TO_SHOW),fontsize='x-small')
	plt.subplot(2,2,3)
	plt.xlabel('Time (min)')
	plt.ylabel('Log10 Absolute Relative Difference')
	plt.subplot(2,2,4)
	plt.xlabel('Time (min)')
	plt.legend(fontsize='xx-small', loc='best')

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
