"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/30/2016
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
import scipy.cluster.hierarchy as sch

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import CMAP_COLORS_255
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells in each seed
		ap = AnalysisPaths(variantDir, cohort_plot = True)

		fig, axesList = plt.subplots(ap.n_seed, ap.n_generation, sharex = False, figsize=(6 + 2*ap.n_generation,10 + 2*ap.n_seed))

		plt.suptitle("Full Metabolic Network Reaction Fluxes")

		for seedNum in xrange(ap.n_seed):
			for generationNum in xrange(ap.n_generation):

				# Only plot one cell per seed per generation
				simDir = ap.get_cells(seed=[seedNum], generation=[generationNum])[0]
				simOutDir = os.path.join(simDir, "simOut")

				fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
				time = fbaResults.readColumn("time")
				reactionFluxes = fbaResults.readColumn("reactionFluxes")
				fbaResults.close()


				normalized = (
					reactionFluxes
					/ (np.mean(np.abs(reactionFluxes), 0) + 2 * np.std(np.abs(reactionFluxes), 0))
					).transpose()

				# Cluster once, and then plot all other cells reactions in that same order
				if seedNum == 0 and generationNum == 0:
					linkage = sch.linkage(reactionFluxes.T)
					linkage[:, 2] = np.fmax(linkage[:, 2], 0) # fixes rounding issues leading to negative distances
					index = sch.leaves_list(linkage)

				cmap = colors.LinearSegmentedColormap.from_list(
					"white to blue with upper extreme",
					CMAP_COLORS
					)

				cmap.set_over(CMAP_OVER)

				norm = colors.Normalize(vmin = -1, vmax = +1)

				if ap.n_seed > 1 and ap.n_generation > 1:
					currentAxes = axesList[seedNum][generationNum]
				elif ap.n_seed > 1 and ap.n_generation == 1:
					currentAxes = axesList[seedNum]
				elif ap.n_seed == 1 and ap.n_generation > 1:
					currentAxes = axesList[generationNum]
				elif ap.n_seed == 1 and ap.n_generation == 1:
					currentAxes = axesList

				currentAxes.imshow(
					normalized[index, :],
					aspect = "auto",
					interpolation='nearest',
					origin = "lower",
					cmap = cmap,
					norm = norm
					)

				means = np.mean(reactionFluxes[:,index],axis=0)

				if seedNum == 0 and generationNum == 0:
					originalMeans = means.copy()

				correlationWithStart = np.corrcoef(originalMeans, means)

				if seedNum == 0 and generationNum == 0:
					currentAxes.set_title("Correlation with seed 0, gen 0: {:.3}".format(correlationWithStart[0,1]), fontsize="x-small")
				else:
					currentAxes.set_title("{:.3}".format(correlationWithStart[0,1]), fontsize="x-small")


				xticks = np.linspace(0,time.size-1, 5, dtype=np.int)
				currentAxes.set_xticks(xticks)

				if generationNum == 0:
					currentAxes.set_ylabel("Seed {}".format(seedNum))

				if seedNum == ap.n_seed - 1:
					if generationNum == 0:
						currentAxes.set_xlabel("Time (min)")
					else:
						currentAxes.set_xlabel("")
					currentAxes.set_xticklabels(np.round(time[xticks]/60.).astype(int))
				else:
					currentAxes.set_xticklabels([])

				currentAxes.set_yticks([])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
