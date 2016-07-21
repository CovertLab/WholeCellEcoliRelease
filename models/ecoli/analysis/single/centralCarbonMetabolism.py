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

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from wholecell.analysis.plotting_tools import CMAP_COLORS_255
from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, VOLUME_UNITS, TIME_UNITS

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

BURN_IN_TIMESTEPS = 500

WITH_FLUX_COLOR = 'blue'
NO_FLUX_COLOR = 'grey'
AVERAGE_COLOR = 'red'

_generatedID_reverseReaction = "{} (reverse)"

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

	coefficients = dryMass / cellMass * cellDensity * (timeStepSec * units.s)

	dryMassFracAverage = np.mean(dryMass / cellMass)

	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])
	toya_fluxes_dict = dict(zip(toya_reactions, toya_fluxes))

	toyaVsReactionAve = []
	for toyaReactionID, toyaFlux in toya_fluxes_dict.iteritems():
		fluxTimeCourse = net_flux(toyaReactionID, reactionIDs, reactionFluxes, reverseRxnFormat=_generatedID_reverseReaction)
		fluxAve = np.mean(fluxTimeCourse)
		toyaVsReactionAve.append((fluxAve.asNumber(FLUX_UNITS), toyaFlux.asNumber(FLUX_UNITS)))

	toyaVsReactionAve = FLUX_UNITS * np.array(toyaVsReactionAve)
	correlationCoefficient = np.corrcoef(toyaVsReactionAve[:,0].asNumber(FLUX_UNITS), toyaVsReactionAve[:,1].asNumber(FLUX_UNITS))[0,1]

	fig = plt.figure(figsize = (30, 15))

	plt.suptitle("Central Carbon Metabolism Reaction Flux vs. Toya 2010 Measured Fluxes", fontsize=22)

	for idx, fluxName in enumerate(toya_reactions):
		reactionTimeCourse = net_flux(fluxName, reactionIDs, reactionFluxes, reverseRxnFormat=_generatedID_reverseReaction).asNumber(FLUX_UNITS).squeeze()
		if reactionTimeCourse[BURN_IN_TIMESTEPS:].any():
			line_color = WITH_FLUX_COLOR
		else:
			line_color = NO_FLUX_COLOR

		ax = plt.subplot(8,4,idx+1)
		ax.plot(time / 60., reactionTimeCourse, linewidth=2, label=fluxName[:32], color=line_color)
		plt.axhline(y=toya_fluxes.asNumber(FLUX_UNITS)[idx], color=AVERAGE_COLOR)
		plt.legend(loc=0)
		plt.xlabel("Time (min)")
		plt.ylabel("Flux ({})".format(FLUX_UNITS.strUnit()))

	ax = plt.subplot(8,4, idx+2)
	ax.plot(0, 0, linewidth=2, label="Zero flux after initial {} time steps".format(BURN_IN_TIMESTEPS), color=NO_FLUX_COLOR)
	ax.plot(0, 0, linewidth=2, label="Nonzero flux".format(BURN_IN_TIMESTEPS), color=WITH_FLUX_COLOR)
	ax.plot(0, 0, linewidth=1, label="Toya 2010 observed flux", color=AVERAGE_COLOR)
	ax.legend(loc = 10)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.set_xticks([])
	ax.set_yticks([])
	# ax.set_title("Grey line denotes zero flux after initial {} time steps.".format(BURN_IN_TIMESTEPS), fontsize=12, bbox={'facecolor':'red', 'alpha':0.5, 'pad':1})

	plt.subplots_adjust()

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

def net_flux(reactionID, reactionIDs, reactionFluxes, reverseRxnFormat="{} (reverse)"):
	fluxTimeCourse = reactionFluxes[:,np.where(reactionIDs == reactionID)]

	reverseID = reverseRxnFormat.format(reactionID)
	if reverseID in reactionIDs:
		reverseTimeCourse = reactionFluxes[:,np.where(reactionIDs == reverseID)]
		fluxTimeCourse = fluxTimeCourse - reverseTimeCourse

	return fluxTimeCourse


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
	
