"""
Plot RNA polymerase counts and counts of mRNA precursors

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
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")

		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnapId = "APORNAP-CPLX[c]"
		rnapIndex = moleculeIds.index(rnapId)
		rnapCountsBulk = bulkMoleculeCounts[:, rnapIndex]

		RNAP_RNA_IDS = ["EG10893_RNA[c]", "EG10894_RNA[c]", "EG10895_RNA[c]", "EG10896_RNA[c]"]
		rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNAP_RNA_IDS], np.int)
		rnapRnaCounts = bulkMoleculeCounts[:, rnapRnaIndexes]

		bulkMolecules.close()

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		rnapIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		nActive = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnapIndex]

		uniqueMoleculeCounts.close()

		plt.figure(figsize = (8.5, 11))

		plt.subplot(5, 1, 1)

		plt.plot(time / 60., nActive + rnapCountsBulk)
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("RNA Polymerase")

		for subplotIdx in xrange(2, 6):
			rnapRnaCountsIdx = subplotIdx - 2

			plt.subplot(5, 1, subplotIdx)

			plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA counts")
			plt.title(RNAP_RNA_IDS[rnapRnaCountsIdx])

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
