#!/usr/bin/env python
"""
Central carbon metabolism comparison to Toya et al for figure 3c

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/13/17
"""

import argparse
import os
import cPickle
import re

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from models.ecoli.analysis.single.centralCarbonMetabolism import net_flux, _generatedID_reverseReaction

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()
	# allDir = ap.get_cells(generation = [0, 1, 2])

	sim_data = cPickle.load(open(simDataFile, "rb"))
	metaboliteNames = np.array(sorted(sim_data.process.metabolism.concDict.keys()))
	nMetabolites = len(metaboliteNames)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
	toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
	toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
	toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

	sim_data = cPickle.load(open(simDataFile))
	cellDensity = sim_data.constants.cellDensity

	modelFluxes = {}
	toyaOrder = []
	for rxn in toyaReactions:
		modelFluxes[rxn] = []
		toyaOrder.append(rxn)

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		mainListener = TableReader(os.path.join(simOutDir, "Main"))
		timeStepSec = mainListener.readColumn("timeStepSec")
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
	for rxn, toyaFlux in toyaFluxesDict.iteritems():
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
	ax = plt.axes()
	ax.set_ylim(plt.xlim())
	whitePadSparklineAxis(plt.axes())

	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
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

	args = parser.parse_args()._Dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
