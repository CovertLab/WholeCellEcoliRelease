"""
Plot NTP counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
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

		NTP_IDS = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
		ntpIndexes = np.array([moleculeIds.index(ntpId) for ntpId in NTP_IDS], np.int)

		ntpCounts = bulkMolecules.readColumn("counts")[:, ntpIndexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 11))

		for idx in xrange(4):

			plt.subplot(2, 2, idx + 1)

			plt.plot(time / 60., ntpCounts[:, idx], linewidth = 2)
			plt.xlabel("Time (min)")
			plt.ylabel("Counts")
			plt.title(NTP_IDS[idx])

		print "NTPs required for cell division (nt/cell-cycle) = %d" % sum(ntpCounts[0, :])
		plt.subplots_adjust(hspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
