"""
Plots cell cycle lengths for all generations.
"""

from __future__ import absolute_import, division, print_function

import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		# Get all cells
		allDir = self.ap.get_cells()

		cellCycleLengths = []
		generations = []
		for idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

			cellCycleLengths.append((time[-1] - time[0]) / 60.)
			generations.append(idx)

		ave_length = np.mean(cellCycleLengths)

		plt.scatter(generations, cellCycleLengths)
		plt.xlabel('Generation')
		plt.ylabel('Time (min)')
		plt.title('Cell cycle lengths\nAverage: {:.1f} min'.format(ave_length))
		plt.axhline(ave_length, color='k', linestyle='--', linewidth=1)
		plt.xticks(generations)
		y_min, y_max = plt.ylim()
		plt.ylim([np.floor(y_min), np.ceil(y_max)])
		plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
