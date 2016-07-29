#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/14/2016
"""

from __future__ import division

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import gridspec

import mpld3
from mpld3 import plugins, utils

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from models.ecoli.analysis.single.centralCarbonMetabolism import net_flux, _generatedID_reverseReaction

from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, VOLUME_UNITS, TIME_UNITS

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
	reactionFluxes = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * np.array(fbaResults.readColumn("reactionFluxes"))
	fbaResults.close()

	fitterPredictedFluxesDict = sim_data.process.metabolism.predictedFluxesDict
	fitterPredictedFluxes = np.array([fitterPredictedFluxesDict[reactionID].asNumber(FLUX_UNITS) for reactionID in reactionIDs])

	# Any seed is fine - this will be consistent within (reruning) a sim, but vary from sim to sim
	np.random.seed(np.abs(hash(simOutDir)) % np.iinfo(np.int32).max)

	n_permutations = 3
	permutedFitterFluxes = np.zeros((n_permutations, fitterPredictedFluxes.shape[0]))
	for iteration in xrange(n_permutations):
		permutation = np.random.permutation(fitterPredictedFluxes)
		permutedFitterFluxes[iteration,:] = permutation

	correlationTimecourse = []
	permutedCorrelationTimecourses = np.zeros((n_permutations, reactionFluxes.asNumber(FLUX_UNITS).shape[0]))
	for timeStep, fluxValue in enumerate(reactionFluxes.asNumber(FLUX_UNITS)):
		correlationCoefficient = np.corrcoef(fluxValue, fitterPredictedFluxes)
		correlationTimecourse.append(correlationCoefficient[0,1])
		for idx, permutedFlux in enumerate(permutedFitterFluxes):
			correlationCoefficient = np.corrcoef(fluxValue, permutedFlux)
			permutedCorrelationTimecourses[idx, timeStep] =  correlationCoefficient[0,1]
	correlationTimecourse = np.array(correlationTimecourse)
	permutedCorrelationTimecourses = np.array(permutedCorrelationTimecourses)

	meanCorrelation = np.nanmean(correlationTimecourse, axis=0)
	permutationMean = np.nanmean(permutedCorrelationTimecourses)

	fig = plt.figure()
	plt.plot(time / 60., correlationTimecourse.T, label="Predicted Fluxes")
	for iteration in xrange(n_permutations):
		plt.plot(time / 60., permutedCorrelationTimecourses[iteration,:], label="Randomly permuted predictions {}".format(iteration + 1))

	plt.title("Fitter Predicted vs. Actual Reaction Fluxes")
	plt.xlabel("Time (min)")
	plt.ylabel("Pearson R")
	plt.legend(loc="best",fontsize='small')
	
	plt.text(.5*np.max(time / 60.),.95*meanCorrelation, "Mean of actual predictions = {:.2}".format(meanCorrelation), horizontalalignment="center")
	plt.text(.5*np.max(time / 60.),np.amax([1.05*np.abs(permutationMean),.03]), "Mean of randomly permuted predictions = {:.2}".format(permutationMean), horizontalalignment="center")

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
	
