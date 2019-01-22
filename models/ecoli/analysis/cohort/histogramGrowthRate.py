from __future__ import division, absolute_import

from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		filepath.makedirs(plotOutDir)

		# Get all cells in each seed
		ap = AnalysisPaths(variantDir, cohort_plot = True)

		max_cells_in_gen = 0
		for genIdx in range(ap.n_generation):
			n_cells = len(ap.get_cells(generation = [genIdx]))

			if n_cells > max_cells_in_gen:
				max_cells_in_gen = n_cells

		fig, axesList = plt.subplots(ap.n_generation,
			sharex=True, sharey=True, figsize=(6, 3*ap.n_generation))

		all_growth_rates = []

		for genIdx in range(ap.n_generation):
			gen_cells = ap.get_cells(generation = [genIdx])
			gen_growth_rates = []

			for simDir in gen_cells:
				try:
					simOutDir = os.path.join(simDir, "simOut")

					# Listeners used
					main_reader = TableReader(os.path.join(simOutDir, "Main"))

					# Load data
					time = main_reader.readColumn("time")
					initialTime = main_reader.readAttribute("initialTime")

					# Get growth rate
					doubling_time = (time[-1] - initialTime) * units.s
					averageGrowthRate = np.log(2) / doubling_time.asNumber(units.min)
					gen_growth_rates.append(averageGrowthRate)
				except Exception:
					# Skip sims that were not able to complete division
					continue

			all_growth_rates.append(np.array(gen_growth_rates))

		# Plot histograms of growth rates
		if ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			gen_growth_rates = all_growth_rates[idx]

			if max_cells_in_gen > 1:
				axes.hist(
					gen_growth_rates,
					int(np.ceil(np.sqrt(len(gen_growth_rates))))
					)
			else:
				axes.plot(gen_growth_rates, 1, 'x')
				axes.set_ylim([0, 2])

			axes.axvline(gen_growth_rates.mean(), color='k', linestyle='dashed', linewidth=2)
			axes.text(
				gen_growth_rates.mean(), 1,
				"Mean = %.4f, Std = %.4f, n = %d" %
				(gen_growth_rates.mean(), gen_growth_rates.std(), len(gen_growth_rates), ))
			axes.set_ylabel("Generation %d" % (idx, ))

		axesList[-1].set_xlabel("Growth rate (1/min)")

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
