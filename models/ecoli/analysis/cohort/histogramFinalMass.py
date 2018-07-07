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

		fig, axesList = plt.subplots(ap.n_generation, sharey = True, sharex = True, subplot_kw=dict((("aspect",0.4),("adjustable",'box-forced'))))

		initial_masses = np.zeros((max_cells_in_gen, ap.n_generation))
		final_masses = np.zeros((max_cells_in_gen, ap.n_generation))

		for genIdx in range(ap.n_generation):
			gen_cells = ap.get_cells(generation = [genIdx])
			for simDir in gen_cells:
				simOutDir = os.path.join(simDir, "simOut")
				mass = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = mass.readColumn("cellMass")

				initial_masses[np.where(simDir == gen_cells)[0], genIdx] = cellMass[0] / 1000.
				final_masses[np.where(simDir == gen_cells)[0], genIdx] = cellMass[-1] / 1000.

		# Plot initial vs final masses
		if ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			if max_cells_in_gen > 1:
				axes.hist(final_masses[:,idx].flatten(), int(np.ceil(np.sqrt(final_masses[:,idx].size))))
			else:
				axes.plot(final_masses[:,idx], 1, 'x')
				axes.set_ylim([0, 2])
			axes.axvline(final_masses[:,idx].mean(), color='k', linestyle='dashed', linewidth=2)
			axes.text(final_masses[:,idx].mean(), 1, "Mean: %.3f Var: %.3f"%(final_masses[:,idx].mean(),final_masses[:,idx].var()))

		axesList[-1].set_xlabel("Initial mass (pg)")
		axesList[ap.n_generation / 2].set_ylabel("Frequency")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
