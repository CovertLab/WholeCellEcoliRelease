#!/usr/bin/env python
"""
Plot reaction max rate over course of the simulation.

@date: Created 7/02/2015
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
import matplotlib.animation as animation

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	
	rateEstimatesArray = enzymeKineticsdata.readColumn("reactionConstraints")
	overconstraintMultiples = enzymeKineticsdata.readColumn("overconstraintMultiples")

	reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()


	fbaData = TableReader(os.path.join(simOutDir, "FBAResults"))

	unconstrainedFluxes = fbaData.readColumn('reactionFluxes')
	fluxNames = fbaData.readAttribute('reactionIDs')
	simulationSteps = fbaData.readColumn('simulationStep')
	fbaData.close()

	# Take the rough halfway point as an example point
	testPoint = len(time)//2

	fluxNamesEstimates = np.array(fluxNames)[np.where(rateEstimatesArray[testPoint] < np.inf)]

	fluxesWithEstimates = unconstrainedFluxes[testPoint][np.where(rateEstimatesArray[testPoint] < np.inf)]
	rateEstimates = rateEstimatesArray[testPoint][np.where(rateEstimatesArray[testPoint] < np.inf)]

	amountOverconstrained = fluxesWithEstimates - rateEstimates
	amountOverconstrained[np.where(amountOverconstrained < 0)] = 0
	overconstrainedFluxes = fluxNamesEstimates[np.where(amountOverconstrained > 0)]
	overconstrainedFluxesDict = dict(zip(overconstrainedFluxes, amountOverconstrained[np.where(amountOverconstrained > 0)]))
	
	relativeRates = fluxesWithEstimates / rateEstimates
	# Only plot rates which are overconstrained
	relativeRates[np.where(relativeRates < 1)] = 0

	plt.figure(figsize=(10,15))

	plt.subplot(3,1,1)
	plt.title("Reaction Fluxes and Kinetic Rates at t={} {}".format(time[testPoint], TIME_UNITS.strUnit()))
	plt.bar(xrange(len(fluxesWithEstimates)),fluxesWithEstimates, color="b", alpha=.7, label="Fluxes")
	plt.bar(xrange(len(rateEstimates)), rateEstimates, .5, color="r", alpha=.5, label="Kinetic Rates")
	plt.xlabel("Reaction")
	plt.ylabel("Reaction Rate ({counts_units}/{volume_units}.{time_units})".format(counts_units=COUNTS_UNITS.strUnit(), volume_units=VOLUME_UNITS.strUnit(), time_units=TIME_UNITS.strUnit()))
	plt.legend(framealpha=.5)

	plt.subplot(3,1,2)
	plt.title("Log Normalized Reaction Fluxes and Kinetic Rates at t={} {}".format(time[testPoint], TIME_UNITS.strUnit()))
	plt.bar(xrange(len(fluxesWithEstimates)),np.log10(fluxesWithEstimates + 1), color="b", alpha=.7, label="Fluxes")
	plt.bar(xrange(len(rateEstimates)), np.log10(rateEstimates + 1), .5, color="r", alpha=.5, label="Kinetic Rates")
	plt.xlabel("Reaction")
	plt.ylabel("Reaction Rate ({counts_units}/{volume_units}.{time_units})".format(counts_units=COUNTS_UNITS.strUnit(), volume_units=VOLUME_UNITS.strUnit(), time_units=TIME_UNITS.strUnit()))
	plt.legend(framealpha=.5)

	plt.subplot(3,1,3)
	plt.title("Kinetic Rates Divided By Reaction Fluxes at t={} {}".format(time[testPoint], TIME_UNITS.strUnit()))

	plt.bar(xrange(len(rateEstimates)),fluxesWithEstimates / rateEstimates, color='purple', alpha=.7, label="Fold Difference")
	plt.xlabel("Reaction")
	plt.ylabel("Fold Difference")
	plt.legend(framealpha=.5)

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
