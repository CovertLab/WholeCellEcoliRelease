"""
Central carbon metabolism comparison to Toya et al for figure 3c

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/13/17
"""

from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle
import re

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from models.ecoli.analysis import multigenAnalysisPlot
import six
from six.moves import zip


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()
		# allDir = ap.get_cells(generation = [0, 1, 2])

		sim_data = cPickle.load(open(simDataFile, "rb"))
		metaboliteNames = np.array(sorted(sim_data.process.metabolism.conc_dict.keys()))
		nMetabolites = len(metaboliteNames)

		validation_data = cPickle.load(open(validationDataFile, "rb"))
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		sim_data = cPickle.load(open(simDataFile, 'rb'))
		cellDensity = sim_data.constants.cellDensity

		modelFluxes = {}
		toyaOrder = []
		for rxn in toyaReactions:
			modelFluxes[rxn] = []
			toyaOrder.append(rxn)

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			mainListener = TableReader(os.path.join(simOutDir, "Main"))
			mainListener.close()

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			massListener.close()

			coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
			reactionFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fbaResults.readColumn("reactionFluxes").T / coefficient).T
			fbaResults.close()

			for toyaReaction in toyaReactions:
				fluxTimeCourse = []

				for rxn in reactionIDs:
					if re.findall(toyaReaction, rxn):
						reverse = 1
						if re.findall("(reverse)", rxn):
							reverse = -1

						if len(fluxTimeCourse):
							fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
						else:
							fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

				if len(fluxTimeCourse):
					modelFluxes[toyaReaction].append(np.mean(fluxTimeCourse).asNumber(units.mmol / units.g / units.h))

		toyaVsReactionAve = []
		for rxn, toyaFlux in six.viewitems(toyaFluxesDict):
			if rxn in modelFluxes:
				toyaVsReactionAve.append((np.mean(modelFluxes[rxn]), toyaFlux.asNumber(units.mmol / units.g / units.h), np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(units.mmol / units.g / units.h)))

		toyaVsReactionAve = np.array(toyaVsReactionAve)
		correlationCoefficient = np.corrcoef(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])[0,1]

		plt.figure(figsize = (8, 8))
		plt.title("Central Carbon Metabolism Flux, Pearson R = {:.2}".format(correlationCoefficient))
		plt.errorbar(toyaVsReactionAve[:,1], toyaVsReactionAve[:,0], xerr = toyaVsReactionAve[:,3], yerr = toyaVsReactionAve[:,2], fmt = "o", ecolor = "k")
		ylim = plt.ylim()
		plt.plot([ylim[0], ylim[1]], [ylim[0], ylim[1]], color = "k")
		plt.xlabel("Toya 2010 Reaction Flux [mmol/g/hr]")
		plt.ylabel("Mean WCM Reaction Flux [mmol/g/hr]")
		ax = plt.gca()
		ax.set_ylim(plt.xlim())
		whitePadSparklineAxis(ax)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
