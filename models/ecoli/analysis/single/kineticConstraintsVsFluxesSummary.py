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

VERBOSE = False

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


	# Only look at fluxes which had a noninfinite estimate at at least one point
	rateEstimates =  rateEstimatesArray[1:]
	rateEstimates += 1.0
	rateEstimates[np.where(rateEstimates == np.inf)] = 0

	fluxesWithEstimates = unconstrainedFluxes[:,np.where(np.sum(rateEstimates, axis=0) > 0)[0]]
	fluxNamesEstimates = np.array([fluxNames[x] for x in np.where(np.sum(rateEstimates, axis=0) > 0)[0]])

	rateEstimates = rateEstimates[:,np.where(np.sum(rateEstimates, axis=0) > 0)[0]]
	rateEstimates -= 1.0
	rateEstimates = np.vstack((np.zeros(rateEstimates.shape[1]),rateEstimates))

	relativeRates = fluxesWithEstimates / rateEstimates

	# relativeRates[np.where(relativeRates < 1)] = 0
	relativeRates[np.isnan(relativeRates)] = 0
	
	means = np.mean(relativeRates, axis=0)
	medians = np.median(relativeRates, axis=0)
	variances = np.var(relativeRates, axis=0)


	overconstraintDict = dict(zip(fluxNamesEstimates, np.mean(relativeRates, axis=0)))
	means = np.mean(relativeRates, axis=0)
	medians = np.median(relativeRates, axis=0)
	variances = np.var(relativeRates, axis=0)

	overconstraintMediansDict = dict(zip(fluxNamesEstimates[np.where(medians > 0)], medians[np.where(medians > 0)]))
	overconstraintMeansDict = dict(zip(fluxNamesEstimates[np.where(means > 0)], means[np.where(means > 0)]))

	if VERBOSE:
		for reaction, overconstraintMultipleAverage in overconstraintDict.iteritems():
			if overconstraintMultipleAverage > 1:
				print reaction, overconstraintMultipleAverage

	plt.figure(figsize=(10,15))

	plt.subplot(2,1,1)
	plt.title("Amount of Overconstraint")
	plt.bar(xrange(relativeRates.shape[1]), means, color="brown", alpha=.7, label="Mean")
	plt.bar(xrange(relativeRates.shape[1]), medians, .5, color="gray", alpha=.8, label="Median")
	plt.xlabel("Reaction")
	plt.ylabel("Rate Estimate / Reaction Flux")
	plt.legend(framealpha=.5)

	plt.subplot(2,1,2)
	plt.title("Log Normalized Amount of Overconstraint")
	plt.bar(xrange(relativeRates.shape[1]),np.log10(means + 1), color="brown", alpha=.7, label="Log10 Mean")
	plt.bar(xrange(relativeRates.shape[1]), np.log10(medians + 1), .5, color="gray", alpha=.8, label=" Log10 Median")
	plt.xlabel("Reaction")
	plt.ylabel("Normalized Log10 (Rate Estimate / Reaction Flux)")
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
