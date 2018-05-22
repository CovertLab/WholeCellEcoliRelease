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
from matplotlib import pyplot as plt
import itertools

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE

BURN_IN_PERIOD = 150

NUMERICAL_ZERO = 1e-15

RANGE_THRESHOLD = 2
MOVING_AVE_WINDOW_SIZE = 200

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

	# Clip reaction fluxes which are less than numerical zero to numerical zero
	reactionFluxes[np.abs(reactionFluxes) < NUMERICAL_ZERO] = 0

	# Build a mapping from reaction to color
	idToColor = {}
	for reactionID, color in itertools.izip(reactionIDs, itertools.cycle(COLORS_LARGE)):
		idToColor[reactionID] = color

	plt.figure(figsize = (17, 11))

	for idx, (reactionID, reactionFlux) in enumerate(zip(reactionIDs, reactionFluxes.T)):
		runningMeanFlux = np.convolve(reactionFlux[BURN_IN_PERIOD:], np.ones((MOVING_AVE_WINDOW_SIZE,))/MOVING_AVE_WINDOW_SIZE, mode='valid')

		meanNormFlux = runningMeanFlux / np.mean(runningMeanFlux)
		fluxRange = meanNormFlux.max() - meanNormFlux.min()

		# Unadjusted
		plt.subplot(2,2,1)
		plt.plot(time / 60., reactionFlux, label=reactionID, color=idToColor[reactionID])

		# Log scale
		plt.subplot(2,2,3)
		plt.plot(time / 60., np.log10(reactionFlux), label=reactionID, color=idToColor[reactionID])


		if fluxRange > RANGE_THRESHOLD:
			# Unadjusted
			plt.subplot(2,2,2)
			plt.plot(time / 60., reactionFlux, label=reactionID, color=idToColor[reactionID])

			# Log scale
			plt.subplot(2,2,4)
			plt.plot(time / 60., np.log10(reactionFlux), label=reactionID, color=idToColor[reactionID])


	plt.suptitle("Reaction Fluxes")
	plt.subplot(2,2,1)
	plt.ylabel('Flux {}'.format(FLUX_UNITS.strUnit()))
	plt.subplot(2,2,2)
	plt.title("Only displaying fluxes whose moving-average (window size {}), range spans at least {}x its mean.".format(MOVING_AVE_WINDOW_SIZE, RANGE_THRESHOLD),fontsize='x-small')
	plt.subplot(2,2,3)
	plt.xlabel('Time (min)')
	plt.ylabel('Log10 Flux {}'.format(FLUX_UNITS.strUnit()))
	plt.subplot(2,2,4)
	plt.xlabel('Time (min)')

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
