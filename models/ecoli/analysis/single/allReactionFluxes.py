"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/14/2016
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt
import itertools

from wholecell.io.tablereader import TableReader

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

BURN_IN_PERIOD = 150

NUMERICAL_ZERO = 1e-15

RANGE_THRESHOLD = 2
MOVING_AVE_WINDOW_SIZE = 200


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

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

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
