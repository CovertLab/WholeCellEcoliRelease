"""
Plots cell cycle lengths for all generations.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/24/2015
"""

from __future__ import absolute_import

import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# Get all cells
		allDir = ap.get_cells()

		cellCycleLengths = []
		generations = []
		for idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

			cellCycleLengths.append((time[-1] - time[0]) / 60.)
			generations.append(idx)

		plt.scatter(generations, cellCycleLengths)
		plt.xlabel('Generation')
		plt.ylabel('Time (min)')
		plt.title('Cell cycle lengths')
		plt.xticks(generations)
		y_min, y_max = plt.ylim()
		plt.ylim([np.floor(y_min), np.ceil(y_max)])
		plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
