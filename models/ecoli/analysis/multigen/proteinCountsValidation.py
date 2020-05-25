"""
Compare protein counts to Schmidt 2015 data set

@author: Javier	Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/4/2017
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, "rb"))
		validation_data = cPickle.load(open(validationDataFile, "rb"))

		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		schmidt_idx = [ids_translation.index(x) for x in validation_data.protein.schmidt2015Data["monomerId"].tolist()]

		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDir = ap.get_cells()

		sim_schmidt_counts_multigen = []

		fig = plt.figure(figsize = (4, 4))

		for simDir in allDir:
			# print(simDir)

			simOutDir = os.path.join(simDir, "simOut")

			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			avgCounts = monomerCounts.readColumn("monomerCounts").mean(axis=0)
			sim_schmidt_counts = avgCounts[schmidt_idx]

			sim_schmidt_counts_multigen.append(sim_schmidt_counts)

		sim_schmidt_counts_multigen = (np.array(sim_schmidt_counts_multigen)).mean(axis = 0)

		axis = plt.subplot(1,1,1)

		axis.plot(np.log10(schmidt_counts + 1), np.log10(sim_schmidt_counts_multigen + 1), 'o', color = "black", markersize = 6, alpha = 0.1, zorder = 1, markeredgewidth = 0.0)
		# print(pearsonr( np.log10(sim_schmidt_counts_mulitgen + 1), np.log10(schmidtCounts + 1) )[0])

		maxLine = np.ceil(
						max((np.log10(schmidt_counts + 1)).max(),
						(np.log10(sim_schmidt_counts_multigen + 1)).max())
					)
		plt.plot([0, maxLine], [0, maxLine], '-k')

		plt.xlim(xmin=0, xmax=maxLine)
		plt.ylim(ymin=0, ymax=maxLine)

		axis.spines["right"].set_visible(False)
		axis.spines["top"].set_visible(False)
		axis.spines["left"].set_position(("outward", 10))
		axis.spines["bottom"].set_position(("outward", 10))
		axis.tick_params(right=False, top=False)
		axis.tick_params(which = "both", direction = "out")

		axis.set_xlim([-0.07, maxLine])
		axis.set_ylim([-0.07, maxLine])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
