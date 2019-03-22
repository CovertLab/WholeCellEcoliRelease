from __future__ import absolute_import

import os

from matplotlib import pyplot as plt

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

		# TODO: Declutter Y-axis

		# Get all cells
		allDir = ap.get_cells().tolist()

		massNames = [
					"dryMass",
					]

		cleanNames = [
					"Dry\nmass",
					]

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			for idx, massType in enumerate(massNames):
				massToPlot = mass.readColumn(massNames[idx])

				f = plt.figure(figsize = (1.25, 0.8), frameon = False)

				ax = f.add_axes([0, 0, 1, 1])
				ax.axis("off")

				ax.plot(time, massToPlot, linewidth = 2)
				ax.set_ylim([massToPlot.min(), massToPlot.max()])
				ax.set_xlim([time.min(), time.max()])

				exportFigure(plt, plotOutDir, "r01_{}_gen{}".format(massType, allDir.index(simDir)))
				plt.close("all")


if __name__ == "__main__":
	Plot().cli()
