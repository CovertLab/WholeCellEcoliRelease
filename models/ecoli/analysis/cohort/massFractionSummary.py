from __future__ import absolute_import

import os
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from wholecell.analysis.plotting_tools import COLORS_256
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

COLORS = [[colorValue/255. for colorValue in color] for color in COLORS_256]


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

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

		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		fig, axesList = plt.subplots(len(massNames), sharex = True)

		currentMaxTime = 0

		# Get all cells in each seed
		ap = AnalysisPaths(variantDir, cohort_plot = True)
		all_cells = ap.get_cells()

		for simDir in all_cells:
			simOutDir = os.path.join(simDir, "simOut")

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			for idx, massType in enumerate(massNames):
				massToPlot = mass.readColumn(massType)
				axesList[idx].plot(((time / 60.) / 60.), massToPlot, linewidth = 2)

				# set axes to size that shows all generations
				cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
				if cellCycleTime > currentMaxTime:
					currentMaxTime = cellCycleTime

				axesList[idx].set_xlim(0, currentMaxTime*int(ap.n_generation)*1.1)
				axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

		for axes in axesList:
			axes.get_ylim()
			axes.set_yticks(list(axes.get_ylim()))

		axesList[0].set_title("Cell mass fractions")
		axesList[len(massNames) - 1].set_xlabel("Time (hr)")
		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
