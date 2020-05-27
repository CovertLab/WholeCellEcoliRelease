"""
Plots the doubling times for all seeds and all generations (after some
threshold number of burn-in generations) as a histogram.

Notes
-----
More strictly speaking, this is the division time or cell-cycle time, for which
are only doubling in an average sense.

There are many hard-coded values in here for the bounds and bins.  This is to
standardize the output across sets of simulations.

"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units

from models.ecoli.analysis import cohortAnalysisPlot

# First generation (counting from zero) from which to gather doubling time
# values.  If fewer generations were run, this script quits early without
# plotting anything.
FIRST_GENERATION = 2

DOUBLING_TIME_BOUNDS_MINUTES = [35, 170]
N_BINS = 30
FREQUENCY_MAX = 160

FIGSIZE = (3.5, 3.5)

HIST_STYLE = dict(
	fc = 'royalblue'
	)

TARGET_LINE_STYLE = dict(
	color = 'crimson',
	lw = 2
	)

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		analysis_paths = AnalysisPaths(variantDir, cohort_plot = True)

		n_gens = analysis_paths.n_generation

		if n_gens - 1 < FIRST_GENERATION:
			print('Not enough generations to plot.')
			return

		sim_dirs = analysis_paths.get_cells(
			generation = range(FIRST_GENERATION, n_gens)
			)

		sim_data = cPickle.load(open(simDataFile, "rb"))

		doubling_times_minutes = []

		for sim_dir in sim_dirs:
			sim_out_dir = os.path.join(sim_dir, "simOut")

			# Assume simulated time == doubling time
			time = TableReader(os.path.join(sim_out_dir, 'Main')).readColumn('time')

			# Time is relative to the first simulation, so need to take a difference
			doubling_time = time[-1] - time[0]

			doubling_times_minutes.append(doubling_time / 60.)

		plt.figure(figsize = FIGSIZE)

		bins = np.linspace(
			DOUBLING_TIME_BOUNDS_MINUTES[0],
			DOUBLING_TIME_BOUNDS_MINUTES[1],
			N_BINS + 1 # +1 because we need n+1 bin bounds for n bins
			)

		plt.hist(doubling_times_minutes, bins = bins, **HIST_STYLE)

		plt.axvline(
			sim_data.doubling_time.asNumber(units.min),
			**TARGET_LINE_STYLE
			)

		plt.title('n = {}'.format(len(doubling_times_minutes)))

		plt.xlim(*DOUBLING_TIME_BOUNDS_MINUTES)
		plt.ylim(0, FREQUENCY_MAX)

		plt.xlabel('Doubling time (minutes)')

		# TODO (John): How to enforce standard axes dimensions?
		# TODO (John): plt.tight_layout()?

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
