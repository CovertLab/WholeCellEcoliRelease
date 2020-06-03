"""
Produces histograms of frequency that production of a metabolite is limited (at least 50 time steps set by WINDOW)

@date: Created 1/12/2017
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

WINDOW = 50


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		fig, axesList = plt.subplots(3)
		fig.set_size_inches(11, 11)

		histo = np.zeros(4)

		ax2 = axesList[2]
		metaboliteNames = None
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			# Listeners used
			enzymeKineticsData = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
			main_reader = TableReader(os.path.join(simOutDir, "Main"))

			# Get names of metabolites and set up to track limited generations
			# in the first sim analyzed
			if metaboliteNames is None:
				metaboliteNames = np.array(enzymeKineticsData.readAttribute("metaboliteNames"))
				nMetabolites = len(metaboliteNames)
				limitedCounts = np.zeros(len(metaboliteNames))

			metaboliteCounts = enzymeKineticsData.readColumn("metaboliteCountsFinal")
			normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]

			# Read time info from the listener
			initialTime = main_reader.readAttribute("initialTime")
			time = main_reader.readColumn("time")

			metaboliteLimited = np.zeros((len(time), nMetabolites))

			diff = np.diff(normalizedCounts, axis = 0)
			limited = []
			for i in xrange(diff.shape[0] - WINDOW):
				currentStepLimited = np.where(np.any(diff[i:i + WINDOW] > 0, axis = 0) == False)[0].astype(int)
				metaboliteLimited[i, currentStepLimited] = 1
				limited = np.append(limited, currentStepLimited).astype(int)

			nLimited = len(np.unique(limited))
			if nLimited >= len(histo):
				histo = np.append(histo, np.zeros(nLimited - len(histo) + 1))
			histo[nLimited] += 1
			limitedCounts[limited] += 1

			ax2.plot(time / 60, metaboliteLimited * range(metaboliteLimited.shape[1]))
			ax2.axvline(initialTime / 60, color = "r", linestyle = "--")

		ax2.set_xlim([0, max(time) / 60])
		ax2.set_xlabel("Time (min)")
		ax2.set_ylabel("Limited")

		ax0 = axesList[0]
		labels = np.arange(len(histo))
		ax0.bar(labels, histo, align='center')
		ax0.set_xticks(labels)
		ax0.set_xlabel("Number of limited metabolites")
		ax0.set_ylabel("Number of generations")

		ax1 = axesList[1]
		ax1.bar(np.arange(len(np.where(limitedCounts > 0)[0])), limitedCounts[limitedCounts > 0], align='center')
		ax1.set_xticks(np.arange(len(np.where(limitedCounts > 0)[0])))
		ax1.set_xticklabels(metaboliteNames[limitedCounts > 0],
			fontsize=6, ha='right', rotation=10)
		ax1.set_xlabel("Metabolite Limited")
		ax1.set_ylabel("Number of genreations")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
