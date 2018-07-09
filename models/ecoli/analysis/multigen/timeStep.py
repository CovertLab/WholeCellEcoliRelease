"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/17/2015
"""

from __future__ import absolute_import
from __future__ import division

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

		# Get all cells
		allDir = ap.get_cells()

		plt.figure(figsize = (8.5, 11))

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			#initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			timeSteps = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			absoluteTime = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			relativeTime = absoluteTime - initialTime

			plt.subplot(3,1,1)
			plt.title("Simulation Time Steps")
			plt.plot(timeSteps)
			plt.xlabel("Increment")
			plt.ylabel("timeStep (s)")

			plt.subplot(3,1,2)
			plt.plot(absoluteTime,timeSteps)
			plt.xlabel("Cell time (s)")
			plt.ylabel("timeStep (s)")

			plt.subplot(3,1,3)
			plt.plot(relativeTime,timeSteps)
			plt.xlabel("Relative cell time (s)")
			plt.ylabel("timeStep (s)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
