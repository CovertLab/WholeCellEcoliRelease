from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells in each seed
		ap = AnalysisPaths(variantDir, cohort_plot = True)

		max_cells_in_gen = 0
		for genIdx in range(ap.n_generation):
			n_cells = len(ap.get_cells(generation = [genIdx]))
			if n_cells > max_cells_in_gen:
				max_cells_in_gen = n_cells

		fig, axesList = plt.subplots(ap.n_generation)

		growth_rate = np.zeros((max_cells_in_gen, ap.n_generation))

		for genIdx in range(ap.n_generation):
			gen_cells = ap.get_cells(generation = [genIdx])
			for simDir in gen_cells:
				simOutDir = os.path.join(simDir, "simOut")
				# growthRate = np.nanmean(TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate"))
				# growthRate = (1 / units.s) * growthRate
				# growthRate = growthRate.asNumber(1 / units.min)

				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
				initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
				doubling_time = (time[-1] - initialTime) * units.s
				averageGrowthRate = np.log(2) / doubling_time.asNumber(units.min)
				growth_rate[np.where(simDir == gen_cells)[0], genIdx] = averageGrowthRate


		norm_growth_rate = growth_rate / growth_rate.mean(axis=0)
		norm_max = norm_growth_rate.max()
		norm_min = norm_growth_rate.min()
		norm_median = np.median(norm_growth_rate)
		norm_std = norm_growth_rate.std()

		bin_max = norm_median + 1*norm_std
		bin_min = norm_median - 1*norm_std
		nbins = 200
		bin_width = (bin_max - bin_min) / nbins

		# Plot initial vs final masses
		if ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			if max_cells_in_gen > 1:
				axes.hist(norm_growth_rate[:,idx].flatten(), bins=np.arange(bin_min, bin_max, bin_width))
			else:
				axes.plot(norm_growth_rate[:,idx], 1, 'x')
				axes.set_ylim([0, 2])

			axes.axvline(norm_growth_rate[:,idx].mean() + norm_growth_rate[:,idx].std(), color='k', linewidth=1, alpha = 0.5)
			axes.axvline(norm_growth_rate[:,idx].mean() - norm_growth_rate[:,idx].std(), color='k', linewidth=1, alpha = 0.5)
			# mean = growth_rate[:,idx].mean()
			# variance = growth_rate[:,idx].var()
			# sd = variance**2.
			# ypos = np.array(axes.get_ylim()).mean()
			#axes.text(growth_rate[:,idx].mean(), ypos, "Mean:{}\nSD:{}\nSD/Mean:{}"%(mean, sd, sd / mean))
			axes.set_xticks([0.5, 0.8, 1.0, 1.2, 1.5])
			axes.set_xlim([bin_min, bin_max])
		axesList[-1].set_xlabel("Normed Growth rate")
		axesList[ap.n_generation / 2].set_ylabel("Frequency")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
