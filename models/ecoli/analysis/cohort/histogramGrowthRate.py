from __future__ import absolute_import, division, print_function

from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader, TableReaderError
from wholecell.utils import units
from six.moves import range


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells in each seed

		max_cells_in_gen = 0
		for genIdx in range(self.ap.n_generation):
			n_cells = len(self.ap.get_cells(generation = [genIdx]))

			if n_cells > max_cells_in_gen:
				max_cells_in_gen = n_cells

		# noinspection PyTypeChecker
		fig, axesList = plt.subplots(self.ap.n_generation,
			sharex=True, sharey=True, figsize=(6, 3*self.ap.n_generation))

		all_growth_rates = []

		for genIdx in range(self.ap.n_generation):
			gen_cells = self.ap.get_cells(generation = [genIdx])
			gen_growth_rates = []

			for simDir in gen_cells:
				simOutDir = os.path.join(simDir, "simOut")
				main_path = os.path.join(simOutDir, "Main")

				try:
					# Listeners used
					main_reader = TableReader(main_path)

					# Load data
					time = main_reader.readColumn("time")
					initialTime = main_reader.readAttribute("initialTime")

					# Get growth rate
					doubling_time = (time[-1] - initialTime) * units.s
					averageGrowthRate = np.log(2) / doubling_time.asNumber(units.min)
					gen_growth_rates.append(averageGrowthRate)
				except (TableReaderError, EnvironmentError) as e:
					# Skip sims that were not able to complete division
					print("Couldn't read the Table {}; maybe the cell didn't finish division; skipping this sim: {!r}"
						.format(main_path, e))
					continue

			all_growth_rates.append(np.array(gen_growth_rates))

		# Plot histograms of growth rates
		if self.ap.n_generation == 1:
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
