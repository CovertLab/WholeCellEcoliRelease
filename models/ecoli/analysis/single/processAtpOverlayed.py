"""
Plot water allocation for each process

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/20/2016
"""

from __future__ import absolute_import
from __future__ import division

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
		processNames = bulkMolecules.readAttribute("processNames")

		atpAllocatedInitial = bulkMolecules.readColumn("atpAllocatedInitial")
		atpRequested = bulkMolecules.readColumn("atpRequested")

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		bulkMolecules.close()


		# Plot
		plt.figure(figsize = (8.5, 11))
		rows = 7
		cols = 2

		for processIndex in np.arange(len(processNames)):
			ax = plt.subplot(rows, cols, processIndex + 1)
			ax.plot(time / 60., atpAllocatedInitial[:, processIndex])
			ax.plot(time / 60., atpRequested[:, processIndex])
			ax.set_title(str(processNames[processIndex]), fontsize = 8, y = 0.85)

			ymin = np.amin([atpAllocatedInitial[:, processIndex], atpRequested[:, processIndex]])
			ymax = np.amax([atpAllocatedInitial[:, processIndex], atpRequested[:, processIndex]])
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			# ax.set_xticks([])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		plt.subplots_adjust(hspace = 2.0, wspace = 2.0)

if __name__ == "__main__":
	Plot().cli()
