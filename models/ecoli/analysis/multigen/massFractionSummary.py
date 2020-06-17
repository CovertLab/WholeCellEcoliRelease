from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# Get all cells
		allDir = ap.get_cells()

		massNames = [
					"dryMass",
					"proteinMass",
					#"tRnaMass",
					"rRnaMass",
					'mRnaMass',
					"dnaMass"
					]

		cleanNames = [
					"Dry\nmass",
					"Protein\nmass",
					#"tRNA\nmass",
					"rRNA\nmass",
					"mRNA\nmass",
					"DNA\nmass"
					]

		# noinspection PyTypeChecker
		fig, axesList = plt.subplots(len(massNames), sharex = True)

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			for idx, massType in enumerate(massNames):
				massToPlot = mass.readColumn(massNames[idx])
				axesList[idx].plot(time / 60. / 60., massToPlot, linewidth = 2)

				axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

		for axes in axesList:
			axes.get_ylim()
			axes.set_yticks(list(axes.get_ylim()))

		axesList[0].set_title("Cell mass fractions")
		axesList[len(massNames) - 1].set_xlabel("Time (hr)")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
