from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot
from six.moves import range


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# TODO: Declutter Y-axis

		# Get first cell from each generation
		firstCellLineage = []
		for gen_idx in range(ap.n_generation):
			firstCellLineage.append(ap.get_cells(generation = [gen_idx])[0])

		massNames = [
					#"dryMass",
					"proteinMass",
					"tRnaMass",
					"rRnaMass",
					'mRnaMass',
					"dnaMass"
					]

		cleanNames = [
					#"Dry\nmass",
					"Protein\nmass frac.",
					"tRNA\nmass frac.",
					"rRNA\nmass frac.",
					"mRNA\nmass frac.",
					"DNA\nmass frac."
					]

		# noinspection PyTypeChecker
		fig, axesList = plt.subplots(len(massNames), sharex = True)

		for simDir in firstCellLineage:
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			massData = np.zeros((len(massNames),time.size))

			for idx, massType in enumerate(massNames):
				massData[idx,:] = mass.readColumn(massNames[idx])

			massData = massData / massData.sum(axis = 0)

			for idx, massType in enumerate(massNames):
				axesList[idx].plot(time / 60, massData[idx,:])
				axesList[idx].set_ylabel(cleanNames[idx])

		for axes in axesList:
			axes.set_yticks(list(axes.get_ylim()))

		axesList[-1].set_xlabel('Time (min)')

		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
