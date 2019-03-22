"""
Plot water count

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/1/2014
"""

from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		moleculeIds = bulkMolecules.readAttribute("objectNames")

		waterIndex = np.array(moleculeIds.index('WATER[c]'), np.int)

		waterCount = bulkMolecules.readColumn("counts")[:, waterIndex]
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		bulkMolecules.close()

		plt.figure(figsize = (8.5, 11))

		plt.plot(time / 60., waterCount, linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("WATER[c] counts")
		plt.title("Counts of water")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
