"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/2016
"""

from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle
import re

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

BURN_IN_TIMESTEPS = 500

WITH_FLUX_COLOR = 'blue'
NO_FLUX_COLOR = 'grey'
AVERAGE_COLOR = 'green'
AVERAGE_COLOR_OPPOSITE_SIGN = 'red'


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		validation_data = cPickle.load(open(validationDataFile, "rb"))
		sim_data = cPickle.load(open(simDataFile, "rb"))

		cellDensity = sim_data.constants.cellDensity

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass") * units.fg
		dryMass = massListener.readColumn("dryMass") * units.fg
		massListener.close()

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
		reactionFluxes = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * np.array(fbaResults.readColumn("reactionFluxes"))
		fbaResults.close()

		dryMassFracAverage = np.mean(dryMass / cellMass)

		toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toya_fluxes = np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])

		plt.figure(figsize = (30, 15))

		plt.suptitle("Central Carbon Metabolism Reaction Flux vs. Toya 2010 Measured Fluxes", fontsize=22)

		idx = -1
		for idx, fluxName in enumerate(toya_reactions):
			reactionTimeCourse = net_flux(fluxName, reactionIDs, reactionFluxes).asNumber(FLUX_UNITS).squeeze()
			if reactionTimeCourse[BURN_IN_TIMESTEPS:].any():
				line_color = WITH_FLUX_COLOR
			else:
				line_color = NO_FLUX_COLOR

			toya_flux = toya_fluxes[idx]
			if np.mean(reactionTimeCourse) / toya_flux > 0:
				ave_line_color = AVERAGE_COLOR
			else:
				ave_line_color = AVERAGE_COLOR_OPPOSITE_SIGN

			ax = plt.subplot(8,4,idx+1)
			ax.plot(time / 60., reactionTimeCourse, linewidth=2, label=fluxName[:32], color=line_color)
			plt.axhline(y=toya_fluxes[idx], color=ave_line_color)
			plt.legend(loc=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Flux ({})".format(FLUX_UNITS.strUnit()))

		ax = plt.subplot(8,4, idx+2)
		ax.plot(0, 0, linewidth=2, label="Zero flux after initial {} time steps".format(BURN_IN_TIMESTEPS), color=NO_FLUX_COLOR)
		ax.plot(0, 0, linewidth=2, label="Nonzero flux", color=WITH_FLUX_COLOR)
		ax.plot(0, 0, linewidth=1, label="Toya 2010 observed flux", color=AVERAGE_COLOR)
		ax.plot(0, 0, linewidth=1, label="Toya 2010 observed flux - opposite sign", color=AVERAGE_COLOR_OPPOSITE_SIGN)
		ax.legend(loc = 10)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.yaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])

		plt.subplots_adjust()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

def net_flux(reactionID, reactionIDs, reactionFluxes):
	fluxTimeCourse = []

	for rxn in reactionIDs:
		if re.findall(reactionID, rxn):
			reverse = 1
			if re.findall("(reverse)", rxn):
				reverse = -1

			if len(fluxTimeCourse):
				fluxTimeCourse += reverse * reactionFluxes[:,
				np.where(reactionIDs == rxn)]
			else:
				fluxTimeCourse = reverse * reactionFluxes[:,
				np.where(reactionIDs == rxn)]

	return fluxTimeCourse


if __name__ == "__main__":
	Plot().cli()
