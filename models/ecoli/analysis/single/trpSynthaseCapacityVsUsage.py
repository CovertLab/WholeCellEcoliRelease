"""
Plot enzymatic capacity of tryptophan synthase vs amount of tryptophan needed by translation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/22/2017
"""

from __future__ import absolute_import

import os
import cPickle

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

BURN_IN = 10


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile, "r"))
		trpIdx = sim_data.moleculeGroups.aaIDs.index("TRP[c]")

		growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))

		trpRequests = growthLimits.readColumn("aaRequestSize")[BURN_IN:, trpIdx]

		growthLimits.close()

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		moleculeIds = bulkMolecules.readAttribute("objectNames")

		trpSynIdx = moleculeIds.index("TRYPSYN[c]")

		trpSynCounts = bulkMolecules.readColumn("counts")[BURN_IN:, trpSynIdx]

		bulkMolecules.close()

		trpSynKcat = 2**( (37. - 25.) / 10.) * 4.1 # From PMID 6402362 (kcat of 4.1/s measured at 25 C)

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")[BURN_IN:] - initialTime
		timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")[BURN_IN:]


		trpSynMaxCapacity = trpSynKcat * trpSynCounts * timeStep

		plt.figure(figsize = (8.5, 11))

		plt.subplot(3, 1, 1)

		plt.plot(time / 60., trpSynMaxCapacity, linewidth = 2)
		plt.ylabel("Tryptophan Synthase Max Capacity")

		plt.subplot(3, 1, 2)

		plt.plot(time / 60., trpRequests, linewidth = 2)
		plt.ylabel("TRP requested by translation")

		plt.subplot(3, 1, 3)

		plt.plot(time / 60., trpSynMaxCapacity / trpRequests, linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("(Max capacity) / (Request)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
