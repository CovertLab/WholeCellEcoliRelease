#!/usr/bin/env python
"""
Compare fluxes in simulation to target fluxes for figure 3C

@date: Created 4/3/17
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os
import cPickle
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.plotting_tools import COLORS_LARGE

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

# ignore data from metabolism burnin period
BURN_IN_TIME = 1

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))

	constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly
	constrainedReactions = np.array(sim_data.process.metabolism.constrainedReactionList)
	useAllConstraints = sim_data.process.metabolism.useAllConstraints
	constraintsToDisable = sim_data.process.metabolism.constraintsToDisable

	targetFluxList = []
	actualFluxList = []
	reactionConstraintlist = []

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		mainListener = TableReader(os.path.join(simOutDir, "Main"))
		time = mainListener.readColumn("time")
		mainListener.close()
		burnIn = time > BURN_IN_TIME
		burnIn[0] = False

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass")
		dryMass = massListener.readColumn("dryMass")
		massListener.close()

		coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		# read constraint data
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		targetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
		actualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
		reactionConstraint = enzymeKineticsReader.readColumn("reactionConstraint")
		enzymeKineticsReader.close()

		targetFluxes = targetFluxes.asNumber(units.mmol / units.g / units.h)
		actualFluxes = actualFluxes.asNumber(units.mmol / units.g / units.h)

		targetAve = np.nanmean(targetFluxes[burnIn, :], axis = 0)
		actualAve = np.nanmean(actualFluxes[burnIn, :], axis = 0)

		if len(targetFluxList) == 0:
			targetFluxList = np.array([targetAve])
			actualFluxList = np.array([actualAve])
			reactionConstraintList = np.array(reactionConstraint[burnIn, :])
		else:
			targetFluxList = np.concatenate((targetFluxList, np.array([targetAve])), axis = 0)
			actualFluxList = np.concatenate((actualFluxList, np.array([actualAve])), axis = 0)
			reactionConstraintList = np.concatenate((reactionConstraintList, np.array(reactionConstraint[burnIn, :])), axis = 0)

	# determine average across all cells
	targetAve = np.nanmean(targetFluxList, axis = 0)
	actualAve = np.nanmean(actualFluxList, axis = 0)

	# categorize reactions that use constraints with only kcat, Km and kcat, or switch between both types of constraints
	kcatOnlyReactions = np.all(constraintIsKcatOnly[reactionConstraintList], axis = 0)
	kmAndKcatReactions = ~np.any(constraintIsKcatOnly[reactionConstraintList], axis = 0)
	mixedReactions = ~(kcatOnlyReactions ^ kmAndKcatReactions)

	disabledReactions = np.zeros(len(constrainedReactions), dtype = bool)
	if not useAllConstraints:
		for rxn in constraintsToDisable:
			idx = np.where(constrainedReactions == rxn)[0]
			disabledReactions[idx] = True
			kcatOnlyReactions[idx] = False
			kmAndKcatReactions[idx] = False
			mixedReactions[idx] = False

	# categorize how well the actual flux matches the target flux
	thresholds = [2, 10]
	categorization = np.zeros(reactionConstraint.shape[1])
	for i, threshold in enumerate(thresholds):
		categorization[kmAndKcatReactions & (actualAve / targetAve < 1. / threshold)] = i + 1
		categorization[actualAve / targetAve > threshold] = i + 1
	categorization[actualAve == 0] = -2
	categorization[actualAve == targetAve] = -1

	# write data for each reaction to a file
	csvFile = open(os.path.join(plotOutDir, plotOutFileName + ".tsv"), "wb")
	output = csv.writer(csvFile, delimiter = "\t")
	output.writerow(["Km and kcat", "Target", "Actual", "Category"])
	for reaction, target, flux, category in zip(constrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], categorization[kmAndKcatReactions]):
		output.writerow([reaction, target, flux, category])

	output.writerow(["kcat only"])
	for reaction, target, flux, category in zip(constrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], categorization[kcatOnlyReactions]):
		output.writerow([reaction, target, flux, category])

	if np.sum(mixedReactions):
		output.writerow(["mixed constraints"])
		for reaction, target, flux, category in zip(constrainedReactions[mixedReactions], targetAve[mixedReactions], actualAve[mixedReactions], categorization[mixedReactions]):
			output.writerow([reaction, target, flux, category])

	if np.sum(disabledReactions):
		output.writerow(["disabled constraints"])
		for reaction, target, flux, category in zip(constrainedReactions[disabledReactions], targetAve[disabledReactions], actualAve[disabledReactions], categorization[disabledReactions]):
			output.writerow([reaction, target, flux, category])

	csvFile.close()

	# add small number to allow plotting of 0 flux on log scale
	targetAve += 1e-6
	actualAve += 1e-6

	# plot data
	plt.figure(figsize = (8, 8))
	ax = plt.axes()
	from scipy.stats import pearsonr
	plt.plot([-8, 4], [-8, 4], 'k')
	plt.plot(np.log10(targetAve[~disabledReactions]), np.log10(actualAve[~disabledReactions]), "ob", markeredgewidth = 0.25, alpha = 0.25)
	plt.xlabel("Log10(Target Flux [mmol/g/hr])")
	plt.ylabel("Log10(Actual Flux [mmol/g/hr])")
	plt.minorticks_off()
	whitePadSparklineAxis(ax)
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	ax.set_yticks(range(int(ylim[0]), int(ylim[1]) + 1, 2))
	ax.set_xticks(range(int(xlim[0]), int(xlim[1]) + 1, 2))

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

	ax.set_xlabel("")
	ax.set_ylabel("")
	ax.set_title("")
	ax.set_xticklabels([])
	ax.set_yticklabels([])

	exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
