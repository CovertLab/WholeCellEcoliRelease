from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot
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
		fig, axesList = plt.subplots(self.ap.n_generation, sharey = True, sharex = True,
			subplot_kw={'aspect': 0.4, 'adjustable': 'box'})

		initial_masses = np.zeros((max_cells_in_gen, self.ap.n_generation))
		final_masses = np.zeros((max_cells_in_gen, self.ap.n_generation))

		for genIdx in range(self.ap.n_generation):
			gen_cells = self.ap.get_cells(generation = [genIdx])
			for simDir in gen_cells:
				simOutDir = os.path.join(simDir, "simOut")
				mass = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = mass.readColumn("cellMass")

				initial_masses[np.where(simDir == gen_cells)[0], genIdx] = cellMass[0] / 1000.
				final_masses[np.where(simDir == gen_cells)[0], genIdx] = cellMass[-1] / 1000.

		# Plot initial vs final masses
		if self.ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			axes.plot(initial_masses[:, idx], final_masses[:, idx], 'o')
			z = np.polyfit(initial_masses[:, idx], final_masses[:, idx], 1)
			p = np.poly1d(z)
			axes.plot(initial_masses[:, idx], p(initial_masses[:, idx]), '--')
			text_x = np.mean(axes.get_xlim())
			text_y = np.mean(axes.get_ylim()) + np.mean(axes.get_ylim())*0.1
			axes.text(text_x, text_y, r"$m_f$=%.3f$\times$$m_i$ + %.3f" % (z[0], z[1]))


		axesList[-1].set_xlabel("Initial mass (pg)")
		axesList[self.ap.n_generation // 2].set_ylabel("Final mass (pg)")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
