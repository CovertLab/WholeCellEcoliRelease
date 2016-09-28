#!/usr/bin/env python
"""
Plot objective reaction components over time

@date: Created 9/27/2016
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import itertools
import warnings

from wholecell.analysis.plotting_tools import COLORS_LARGE

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

BURN_IN_PERIOD = 150
MAX_NUM_DEVIANTS_TO_SHOW = 5

NUMERICAL_ZERO = 1e-20

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaData = TableReader(os.path.join(simOutDir, "FBAResults"))
	
	kineticTargetFluxNames = fbaData.readAttribute("kineticTargetFluxNames")
	kineticObjectiveComponents = fbaData.readColumn("kineticObjectiveValues")
	homeostaticObjectiveWeights = fbaData.readColumn("homeostaticObjectiveWeight")

	kineticObjectiveWeights = np.ones_like(homeostaticObjectiveWeights) - homeostaticObjectiveWeights

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fbaData.close()


	# Build a mapping from molecule to color
	idToColor = {}
	for fluxID, color in itertools.izip(kineticTargetFluxNames, itertools.cycle(COLORS_LARGE)):
		idToColor[fluxID] = color
	
	plt.figure(figsize = (17, 11))

	# Determine where to set the threshold to hit the target number of most deviant fluxes
	with warnings.catch_warnings() as w:
		warnings.simplefilter("ignore")
		maxima = np.nanmax(kineticObjectiveComponents[BURN_IN_PERIOD:], axis=0)
	maxima = np.delete(maxima, np.where(np.isnan(maxima)))
	threshold_value = sorted(maxima)[-(MAX_NUM_DEVIANTS_TO_SHOW+1)]

	alwaysZero = set()
	for idx, (fluxID, objectiveValues) in enumerate(zip(kineticTargetFluxNames, kineticObjectiveComponents.T)):

		logObjectiveValues = objectiveValues.copy()
		logObjectiveValues[logObjectiveValues == 0] = NUMERICAL_ZERO
		logObjectiveValues = np.log10(logObjectiveValues)

		# Unadjusted
		plt.subplot(2,2,1)
		plt.plot(time / 60., objectiveValues, label=fluxID, color=idToColor[fluxID])

		# Log scale
		plt.subplot(2,2,3)
		plt.plot(time / 60., logObjectiveValues, label=fluxID, color=idToColor[fluxID])

		if ((np.abs(objectiveValues[BURN_IN_PERIOD:]) > threshold_value).any()):
			# Unadjusted
			plt.subplot(2,2,2)
			plt.plot(time / 60., objectiveValues, label=fluxID, color=idToColor[fluxID])

			# Log scale
			plt.subplot(2,2,4)
			plt.plot(time / 60., logObjectiveValues, label=fluxID, color=idToColor[fluxID])

		if not objectiveValues[BURN_IN_PERIOD:].any():
			alwaysZero.add(fluxID)

	plt.suptitle("Metabolism Kinetic Objective Components: {} ({:.1f}%) equal zero after a {} step burn-in.".format(len(alwaysZero), 100*len(alwaysZero)/len(kineticTargetFluxNames), BURN_IN_PERIOD))
	plt.subplot(2,2,1)
	plt.ylabel('Objective value')
	plt.subplot(2,2,2)
	plt.title("Only displaying fluxes with values greater than 0 after a {} step burn-in.".format(BURN_IN_PERIOD),fontsize='x-small')
	plt.title("Top {} fluxes by max objective value after {} step burn-in.".format(MAX_NUM_DEVIANTS_TO_SHOW, BURN_IN_PERIOD),fontsize='x-small')

	plt.subplot(2,2,3)
	plt.xlabel('Time (min)')
	plt.ylabel('Log10 objective value')
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
