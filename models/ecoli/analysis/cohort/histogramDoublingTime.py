from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
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

		fig, axesList = plt.subplots(ap.n_generation, sharex = True)

		doubling_time = np.zeros((max_cells_in_gen, ap.n_generation))

		for genIdx in range(ap.n_generation):
			gen_cells = ap.get_cells(generation = [genIdx])
			for simDir in gen_cells:
				simOutDir = os.path.join(simDir, "simOut")
				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
				initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")

				doubling_time[np.where(simDir == gen_cells)[0], genIdx] = (time.max() - initialTime) / 60.

		# Plot initial vs final masses
		if ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			if max_cells_in_gen > 1:
				axes.hist(doubling_time[:,idx].flatten(), int(np.ceil(np.sqrt(doubling_time[:,idx].size))))
			else:
				axes.plot(doubling_time[:,idx], 1, 'x')
				axes.set_ylim([0, 2])
			axes.axvline(doubling_time[:,idx].mean(), color='k', linestyle='dashed', linewidth=2)
			axes.text(doubling_time[:,idx].mean(), 1, "Mean: %.3f Var: %.3f"%(doubling_time[:,idx].mean(),doubling_time[:,idx].var()))

		axesList[-1].set_xlabel("Doubling time (min))")
		axesList[ap.n_generation / 2].set_ylabel("Frequency")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
