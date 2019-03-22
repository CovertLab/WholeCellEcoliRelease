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

from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units

from models.ecoli.analysis import cohortAnalysisPlot

THROW_ON_BAD_SIMULATION_OUTPUT = False

# First generation (counting from zero) from which to gather doubling time
# values.  If fewer generations were run, this script quits early without
# plotting anything.
FIRST_GENERATION = 2

DOUBLING_TIME_BOUNDS_MINUTES = [0, 170]
N_BINS = 30
FREQUENCY_MAX = 350

FIGSIZE = (3.5, 3.5)

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		analysis_paths = AnalysisPaths(variantDir, cohort_plot = True)

		n_gens = analysis_paths.n_generation

		if n_gens - 1 < FIRST_GENERATION:
			print 'Not enough generations to plot.'
			return

		sim_dirs = analysis_paths.get_cells(
			generation = range(FIRST_GENERATION, n_gens)
			)

		sim_data = cPickle.load(open(simDataFile, "rb"))

		doubling_times_minutes = []

		missing_files = []
		broken_files = []

		for sim_dir in sim_dirs:
			sim_out_dir = os.path.join(sim_dir, "simOut")

			path = os.path.join(sim_out_dir, 'Main')

			if not os.path.exists(path):
				missing_files.append(path)
				continue

			# Assume simulated time == doubling time
			try:
				time = TableReader(path).readColumn('time')

			except Exception as e:
				broken_files.append(path)
				continue

			# Time is relative to the first simulation, so need to take a difference
			try:
				doubling_time = time[-1] - time[0]

			except Exception as e:
				broken_files.append(path)
				continue

			doubling_times_minutes.append(doubling_time / 60.)

		if missing_files or broken_files:
			messages = []

			if missing_files:
				messages.append('Missing files:\n{}'.format(
					'\n'.join(missing_files)
					))

			if broken_files:
				messages.append('Broken files:\n{}'.format(
					'\n'.join(broken_files)
					))

			message = '\n'.join(messages)

			if THROW_ON_BAD_SIMULATION_OUTPUT:
				# Throw late so we get a full picture of what files are missing
				raise Exception(message)

			else:
				print message

		plt.figure(figsize = FIGSIZE)
		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

		bins = np.linspace(
			DOUBLING_TIME_BOUNDS_MINUTES[0],
			DOUBLING_TIME_BOUNDS_MINUTES[1],
			N_BINS + 1 # +1 because we need n+1 bin bounds for n bins
			)

		plt.hist(doubling_times_minutes, bins = bins)

		plt.axvline(
			np.median(doubling_times_minutes), color='k', lw=2, linestyle='--'
			)

		plt.axvline(
			sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min),
			color=color_cycle[2], lw=2
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
