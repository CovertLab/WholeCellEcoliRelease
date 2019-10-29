"""
Plots the doubling times for all seeds and all generations (after some
threshold number of burn-in generations) as a violin plot.

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

DOUBLING_TIME_BOUNDS_MINUTES = [0, 200]
FIGSIZE = (2.5, 5)

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
			generation = range(FIRST_GENERATION, n_gens), seed = range(8))

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

		if not len(doubling_times_minutes):
			print('Skipping plot due to no viable sims.')
			return

		fig, ax = plt.subplots(1, 1, figsize=FIGSIZE)
		ax.violinplot(doubling_times_minutes)
		ax.axhline(
			sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min),
			color='r', lw=1)
		ax.set_ylim(*DOUBLING_TIME_BOUNDS_MINUTES)
		ax.set_xlim([0.5, 1.5])
		ax.set_xticks([])
		y_ticks = ax.get_yticks()
		ax.set_yticklabels([])
		ax.spines['right'].set_visible(False)
		exportFigure(plt, plotOutDir, '{}__clean'.format(plotOutFileName), None)

		ax.set_title('n = {}'.format(len(doubling_times_minutes)))
		ax.set_ylabel('Doubling time (minutes)')
		ax.set_yticks(y_ticks)
		ax.spines['right'].set_visible(True)
		ax.spines['left'].set_visible(True)
		plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
