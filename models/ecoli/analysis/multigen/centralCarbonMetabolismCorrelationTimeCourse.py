#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/2016
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
from scipy.stats import pearsonr

import mpld3
from mpld3 import plugins, utils

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from models.ecoli.analysis.single.centralCarbonMetabolism import net_flux, _generatedID_reverseReaction

from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, VOLUME_UNITS, TIME_UNITS

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)
	
	validation_data = cPickle.load(open(validationDataFile, "rb"))
	sim_data = cPickle.load(open(simDataFile, "rb"))

	cellDensity = sim_data.constants.cellDensity

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# Get all cells
	allDir = ap.get_cells()

	plt.figure(figsize = (8.5, 11))

	currentMaxTime = 0

	fig = plt.figure()
	fig.hold(True)

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
	
		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass") * units.fg
		dryMass = massListener.readColumn("dryMass") * units.fg
		massListener.close()

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
		reactionFluxes = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * np.array(fbaResults.readColumn("reactionFluxes"))
		fluxes_dict = dict(zip(reactionIDs, reactionFluxes))
		fbaResults.close()
		
		coefficients = dryMass / cellMass * cellDensity * (timeStepSec * units.s)
		dryMassFracAverage = np.mean(dryMass / cellMass)

		toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])

		netFluxes = []
		for toyaReactionID in toya_reactions:
			if toyaReactionID in reactionIDs:
				fluxTimeCourse = net_flux(toyaReactionID, reactionIDs, reactionFluxes, reverseRxnFormat=_generatedID_reverseReaction).asNumber(FLUX_UNITS).squeeze()
				netFluxes.append(fluxTimeCourse)

		trimmedReactions = FLUX_UNITS * np.array(netFluxes)

		corrCoefTimecourse = []
		for fluxes in trimmedReactions.asNumber(FLUX_UNITS).T:
			correlationCoefficient = np.corrcoef(fluxes, toya_fluxes.asNumber(FLUX_UNITS))[0,1]
			corrCoefTimecourse.append(correlationCoefficient)

		plt.plot(time / 60., corrCoefTimecourse)
		plt.title("Measured vs. Simulated Central Carbon Fluxes")
		plt.xlabel("Time (min)")
		plt.ylabel("Pearson R")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
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










	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])

	netFluxes = []
	for toyaReactionID in toya_reactions:
		if toyaReactionID in reactionIDs:
			fluxTimeCourse = net_flux(toyaReactionID, reactionIDs, reactionFluxes, reverseRxnFormat=_generatedID_reverseReaction).asNumber(FLUX_UNITS).squeeze()
			netFluxes.append(fluxTimeCourse)

	trimmedReactions = FLUX_UNITS * np.array(netFluxes)

	corrCoefTimecourse = []
	for fluxes in trimmedReactions.asNumber(FLUX_UNITS).T:
		correlationCoefficient = np.corrcoef(fluxes, toya_fluxes.asNumber(FLUX_UNITS))[0,1]
		corrCoefTimecourse.append(correlationCoefficient)

	meanCorr = np.mean(np.array(corrCoefTimecourse)[~np.isnan(corrCoefTimecourse)])

	fig = plt.figure()
	plt.plot(time / 60., corrCoefTimecourse)
	plt.axhline(y=meanCorr, color='r')
	plt.title("Measured vs. Simulated Central Carbon Fluxes")
	plt.text(.5*np.max(time / 60.),.95*meanCorr, "Mean = {:.2}".format(meanCorr), horizontalalignment="center")
	plt.xlabel("Time (min)")
	plt.ylabel("Pearson R")

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
	