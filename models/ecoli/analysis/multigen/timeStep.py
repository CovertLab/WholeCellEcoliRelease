from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		# Get all cells
		allDir = self.ap.get_cells()

		plt.figure(figsize = (8.5, 11))

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			timeSteps = main_reader.readColumn("timeStepSec")
			initialTime = main_reader.readAttribute("initialTime")
			absoluteTime = main_reader.readColumn("time")
			relativeTime = absoluteTime - initialTime

			self.subplot(3,1,1)
			plt.title("Simulation Time Steps")
			plt.plot(timeSteps)
			plt.xlabel("Increment")
			plt.ylabel("timeStep (s)")

			self.subplot(3,1,2)
			plt.plot(absoluteTime,timeSteps)
			plt.xlabel("Cell time (s)")
			plt.ylabel("timeStep (s)")

			self.subplot(3,1,3)
			plt.plot(relativeTime,timeSteps)
			plt.xlabel("Relative cell time (s)")
			plt.ylabel("timeStep (s)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
