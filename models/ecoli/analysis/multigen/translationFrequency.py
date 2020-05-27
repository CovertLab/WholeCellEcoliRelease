"""
Plots frequency of observing at least 1 protein.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/11/2017
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		proteinPresence = []
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			# Get boolean protein presence
			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			proteinCounts = monomerCounts.readColumn("monomerCounts")
			meanProteinCounts = proteinCounts.mean(axis=0)

			proteinPresence.append(meanProteinCounts != 0)

		proteinPresence = np.array(proteinPresence)

		# Plot
		fig = plt.figure(figsize = (12, 12))
		ax = plt.subplot(1, 1, 1)
		nGens = len(allDir)
		ax.hist(np.mean(proteinPresence, axis = 0), nGens)
		ax.set_xlabel("Frequency of observing at least 1 protein copy in 1 generation", fontsize = 14)
		ax.set_ylabel("Number of proteins", fontsize = 14)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
