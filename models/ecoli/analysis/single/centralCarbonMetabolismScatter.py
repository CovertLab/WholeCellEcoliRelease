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
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from scipy.stats import pearsonr

import mpld3
from mpld3 import plugins, utils

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from models.ecoli.analysis.single.centralCarbonMetabolism import net_flux, _generatedID_reverseReaction
from wholecell.analysis.plotting_tools import CMAP_COLORS_255

from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, VOLUME_UNITS, TIME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	sim_data = cPickle.load(open(simDataFile, "rb"))

	cellDensity = sim_data.constants.cellDensity

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
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

	dryMassFracAverage = np.mean(dryMass / cellMass)

	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])
	toya_stdev = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]])
	toya_fluxes_dict = dict(zip(toya_reactions, toya_fluxes))
	toya_stdev_dict = dict(zip(toya_reactions, toya_stdev))

	toyaVsReactionAve = []
	toya_order = []
	for toyaReactionID, toyaFlux in toya_fluxes_dict.iteritems():
		fluxTimeCourse = []

		for rxn in reactionIDs:
			if re.findall(toyaReactionID, rxn):
				reverse = 1
				if re.findall("(reverse)", rxn):
					reverse = -1

				if len(fluxTimeCourse):
					fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
				else:
					fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

		if len(fluxTimeCourse):
			fluxAve = np.mean(fluxTimeCourse)
			fluxStdev = np.std(fluxTimeCourse.asNumber(FLUX_UNITS))
			toyaVsReactionAve.append((fluxAve.asNumber(FLUX_UNITS), toyaFlux.asNumber(FLUX_UNITS), fluxStdev, toya_stdev_dict[toyaReactionID].asNumber(FLUX_UNITS)))
			toya_order.append(toyaReactionID)

	toyaVsReactionAve = FLUX_UNITS * np.array(toyaVsReactionAve)
	correlationCoefficient = np.corrcoef(toyaVsReactionAve[:,0].asNumber(FLUX_UNITS), toyaVsReactionAve[:,1].asNumber(FLUX_UNITS))[0,1]

	fig = plt.figure()

	plt.title("Central Carbon Metabolism Flux, Pearson R = {:.2}".format(correlationCoefficient))
	points = plt.errorbar(toyaVsReactionAve[:,1].asNumber(FLUX_UNITS), toyaVsReactionAve[:,0].asNumber(FLUX_UNITS), xerr = toyaVsReactionAve[:,3].asNumber(FLUX_UNITS),
		yerr = toyaVsReactionAve[:,2].asNumber(FLUX_UNITS), fmt = "o", ecolor = "k")
	plt.ylabel("Mean WCM Reaction Flux {}".format(FLUX_UNITS.strUnit()))
	plt.xlabel("Toya 2010 Reaction Flux {}".format(FLUX_UNITS.strUnit()))

	labels = list(toya_order)
	tooltip = plugins.PointLabelTooltip(points, labels)
	plugins.connect(fig, tooltip)

	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)
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
	
