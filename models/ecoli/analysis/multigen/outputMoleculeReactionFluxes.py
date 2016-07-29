#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2016
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import itertools

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE

BURN_IN_PERIOD = 150
RANGE_THRESHOLD = 1
MOVING_AVE_WINDOW_SIZE = 300

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# Get all cells
	allDir = ap.get_cells()

	plt.figure(figsize = (17, 11))
	idToColor = {}
	for idx, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		outputFluxes = fbaResults.readColumn("outputFluxes")
		outputMoleculeIDs = np.array(fbaResults.readAttribute("outputMoleculeIDs"))
		fbaResults.close()

		# Build a color mapping with the first cell, then stick to it
		if idx == 0:
			idToColor = {}
			for outputMoleculeID, color in itertools.izip(outputMoleculeIDs, itertools.cycle(COLORS_LARGE)):
				idToColor[outputMoleculeID] = color

		for idx, (outputMoleculeID, outputFlux) in enumerate(zip(outputMoleculeIDs, outputFluxes.T)):
			runningMeanFlux = np.convolve(outputFlux[BURN_IN_PERIOD:], np.ones((MOVING_AVE_WINDOW_SIZE,))/MOVING_AVE_WINDOW_SIZE, mode='valid')

			meanNormFlux = runningMeanFlux / np.mean(runningMeanFlux)
			fluxRange = meanNormFlux.max() - meanNormFlux.min()

			# Unadjusted
			plt.subplot(2,2,1)
			plt.plot(time / 60. / 60., outputFlux, label=outputMoleculeID, color=idToColor[outputMoleculeID])

			# Log scale
			plt.subplot(2,2,3)
			plt.plot(time / 60. / 60., np.log10(outputFlux), label=outputMoleculeID, color=idToColor[outputMoleculeID])

			if fluxRange > RANGE_THRESHOLD:
				# Unadjusted
				plt.subplot(2,2,2)
				plt.plot(time / 60. / 60., outputFlux, label=outputMoleculeID, color=idToColor[outputMoleculeID])

				# Log scale
				plt.subplot(2,2,4)
				plt.plot(time / 60. / 60., np.log10(outputFlux), label=outputMoleculeID, color=idToColor[outputMoleculeID])

	plt.suptitle("Output Reaction Fluxes")
	plt.subplot(2,2,1)
	plt.ylabel('Flux {}'.format(FLUX_UNITS.strUnit()))
	plt.subplot(2,2,2)
	plt.title("Only displaying fluxes whose moving-average (window size {}), range spans at least {}x its mean.".format(MOVING_AVE_WINDOW_SIZE, RANGE_THRESHOLD),fontsize='x-small')
	plt.subplot(2,2,3)
	plt.xlabel('Time (min)')
	plt.ylabel('Log10 Flux {}'.format(FLUX_UNITS.strUnit()))
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
