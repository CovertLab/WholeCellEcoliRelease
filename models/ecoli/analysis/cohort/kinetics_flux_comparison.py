"""
Compare fluxes in simulation to target fluxes

@date: Created 4/3/17
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle
import csv

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

# ignore data from metabolism burnin period
BURN_IN_TIME = 1


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(variantDir, cohort_plot = True)
		allDir = ap.get_cells()

		sim_data = cPickle.load(open(simDataFile, "rb"))

		constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly

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
			constrainedReactions = np.array(enzymeKineticsReader.readAttribute("constrainedReactions"))
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

		# categorize how well the actual flux matches the target flux
		thresholds = [2, 10]
		categorization = np.zeros(reactionConstraint.shape[1])
		for i, threshold in enumerate(thresholds):
			categorization[actualAve / targetAve < 1. / threshold] = i + 1
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

		csvFile.close()

		# add small number to allow plotting of 0 flux on log scale
		targetAve += 1e-6
		actualAve += 1e-6

		pearsonAll = pearsonr(np.log10(targetAve), np.log10(actualAve))
		pearsonNoZeros = pearsonr(np.log10(targetAve[(categorization != -2)]), np.log10(actualAve[(categorization != -2)]))

		# plot data
		plt.figure(figsize = (4, 4))
		ax = plt.axes()
		plt.plot([-6, 4], [-6, 4], 'k', linewidth = 0.75)
		plt.plot([-5, 4], [-6, 3], 'k', linewidth = 0.5)
		plt.plot([-6, 3], [-5, 4], 'k', linewidth = 0.5)
		plt.plot(np.log10(targetAve), np.log10(actualAve), 'o', color = "black", markersize = 8, alpha = 0.15, zorder=1, markeredgewidth = 0.0)
		plt.xlabel("Log10(Target Flux [mmol/g/hr])")
		plt.ylabel("Log10(Actual Flux [mmol/g/hr])")
		plt.title("PCC = %.3f, p = %s\n(%.3f, p = %s without points at zero)" % (pearsonAll[0], pearsonAll[1], pearsonNoZeros[0], pearsonNoZeros[1]))
		plt.minorticks_off()
		whitePadSparklineAxis(ax)
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_ylim(ylim[0] - 0.5, ylim[1])
		ax.set_xlim(xlim[0] - 0.5, xlim[1])
		ax.set_yticks(range(-6, int(ylim[1]) + 1, 2))
		ax.set_xticks(range(-6, int(xlim[1]) + 1, 2))

		exportFigure(plt, plotOutDir, plotOutFileName)

		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_title("")
		ax.set_xticklabels([])
		ax.set_yticklabels([])

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
