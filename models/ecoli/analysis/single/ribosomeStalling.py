"""
Plot ribosome stalling

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/22/2014
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

import os

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

		ribosomeData = TableReader(os.path.join(simOutDir, "RibosomeData"))

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		fractionStalled = ribosomeData.readColumn("fractionStalled")

		ribosomeData.close()

		plt.figure(figsize = (8.5, 11))
		plt.plot(time / 60, fractionStalled)

		plt.xlabel("Time (min)")
		plt.ylabel("Fraction of ribosomes stalled")

		plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
