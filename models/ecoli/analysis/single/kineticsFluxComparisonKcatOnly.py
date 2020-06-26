"""
Compare fluxes in simulation to target fluxes

@date: Created 12/15/16
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure


BURN_IN_STEPS = 20


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# read constraint data
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		allTargetFluxes = enzymeKineticsReader.readColumn("targetFluxes")
		allActualFluxes = enzymeKineticsReader.readColumn("actualFluxes")
		kineticsConstrainedReactions = np.array(enzymeKineticsReader.readAttribute("kineticsConstrainedReactions"))
		constraint_is_kcat_only = np.array(enzymeKineticsReader.readAttribute('constraint_is_kcat_only'))

		# kinetic target fluxes
		targetFluxes = allTargetFluxes[:, 0:len(kineticsConstrainedReactions)]
		actualFluxes = allActualFluxes[:, 0:len(kineticsConstrainedReactions)]

		targetAve = np.mean(targetFluxes[BURN_IN_STEPS:, :], axis = 0)
		actualAve = np.mean(actualFluxes[BURN_IN_STEPS:, :], axis = 0)

		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only

		kmAndKcatThresholds = [2, 10]
		kmAndKcatCategorization = np.zeros(np.sum(kmAndKcatReactions))
		for i, threshold in enumerate(kmAndKcatThresholds):
			kmAndKcatCategorization[targetAve[kmAndKcatReactions] / actualAve[kmAndKcatReactions] > threshold] = i + 1
			kmAndKcatCategorization[actualAve[kmAndKcatReactions] / targetAve[kmAndKcatReactions] > threshold] = i + 1
		kmAndKcatCategorization[actualAve[kmAndKcatReactions] == 0] = -1

		kcatOnlyThresholds = [2, 10]
		kcatOnlyCategorization = np.zeros(np.sum(kcatOnlyReactions))
		for i, threshold in enumerate(kcatOnlyThresholds):
			kcatOnlyCategorization[actualAve[kcatOnlyReactions] / targetAve[kcatOnlyReactions] > threshold] = i + 1
		kcatOnlyCategorization[actualAve[kcatOnlyReactions] == 0] = -1

		targetAve += 1e-13
		actualAve += 1e-13

		plt.figure(figsize = (8, 8))
		plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 0], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 0], "og")
		plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 1], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 1], "o")
		plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 2], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 2], "or")
		plt.loglog([1e-12, 1], [1e-12, 1], '--g')
		plt.loglog([1e-12, 1], [1e-11, 10], '--r')
		plt.xlabel("Target Flux (dmol/L/s)")
		plt.ylabel("Actual Flux (dmol/L/s)")

		exportFigure(plt, plotOutDir, plotOutFileName)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
