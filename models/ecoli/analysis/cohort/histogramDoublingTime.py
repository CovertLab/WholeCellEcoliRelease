from __future__ import division, absolute_import

from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		filepath.makedirs(plotOutDir)

		# Get all cells in each seed
		ap = AnalysisPaths(variantDir, cohort_plot=True)

		max_cells_in_gen = 0
		for genIdx in range(ap.n_generation):
			n_cells = len(ap.get_cells(generation = [genIdx]))

			if n_cells > max_cells_in_gen:
				max_cells_in_gen = n_cells

		fig, axesList = plt.subplots(ap.n_generation,
			sharex=True, sharey=True, figsize=(6, 3*ap.n_generation))

		all_doubling_times = []

		for genIdx in range(ap.n_generation):
			gen_cells = ap.get_cells(generation = [genIdx])
			gen_doubling_times = []

			for simDir in gen_cells:
				try:
					simOutDir = os.path.join(simDir, "simOut")

					# Listeners used
					main_reader = TableReader(os.path.join(simOutDir, "Main"))

					# Load data
					time = main_reader.readColumn("time")
					initial_time = main_reader.readAttribute("initialTime")

					# Get doubling time
					gen_doubling_times.append((time[-1] - initial_time) / 60.)
				except Exception:
					# Skip sims that were not able to complete division
					continue

			all_doubling_times.append(np.array(gen_doubling_times))

		# Plot histograms of doubling times
		if ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			gen_doubling_times = all_doubling_times[idx]

			if max_cells_in_gen > 1:
				axes.hist(
					gen_doubling_times,
					int(np.ceil(np.sqrt(len(gen_doubling_times))))
					)
			else:
				axes.plot(gen_doubling_times, 1, 'x')
				axes.set_ylim([0, 2])

			axes.axvline(gen_doubling_times.mean(), color='k', linestyle='dashed', linewidth=2)
			axes.text(
				gen_doubling_times.mean(), 1,
				"Mean = %.2f, Std = %.2f, n = %d" %
				(gen_doubling_times.mean(), gen_doubling_times.std(), len(gen_doubling_times), ))
			axes.set_ylabel("Generation %d" % (idx, ))

		axesList[-1].set_xlabel("Doubling time (min))")

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
