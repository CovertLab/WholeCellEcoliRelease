from __future__ import absolute_import, division, print_function

import os
import itertools

from matplotlib import pyplot as plt
from six.moves import zip

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
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

		ap = AnalysisPaths(inputDir, variant_plot = True)
		all_cells = ap.get_cells()

		# Build a mapping from variant id to color
		idToColor = {}
		for idx, (cell_id, color) in enumerate(zip(all_cells, itertools.cycle(COLORS_LARGE))):
			idToColor[idx] = color

		# noinspection PyTypeChecker
		fig, axesList = plt.subplots(len(massNames), sharex = True)

		currentMaxTime = 0
		for cellIdx, simDir in enumerate(all_cells):
			with open(os.path.join(simDir[:-32],'metadata','short_name')) as f:
				variant_name = [line for line in f][0]

			simOutDir = os.path.join(simDir, "simOut")

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			for massIdx, massType in enumerate(massNames):
				massToPlot = mass.readColumn(massType)
				axesList[massIdx].plot(((time / 60.) / 60.), massToPlot, linewidth = 2, color=idToColor[cellIdx], label=variant_name)

				# set axes to size that shows all generations
				cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
				if cellCycleTime > currentMaxTime:
					currentMaxTime = cellCycleTime

				axesList[massIdx].set_xlim(0, currentMaxTime*ap.n_generation*1.1)
				axesList[massIdx].set_ylabel(cleanNames[massIdx] + " (fg)")

		for idx, axes in enumerate(axesList):
			axes.get_ylim()
			axes.set_yticks(list(axes.get_ylim()))

		axesList[0].set_title("Cell mass fractions")
		plt.legend(bbox_to_anchor=(.92, 5), loc=2, borderaxespad=0., prop={'size':6})
		axesList[len(massNames) - 1].set_xlabel("Time (hr)")
		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
