"""
Plot RNA polymerase counts and counts of mRNA precursors

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import absolute_import

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

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnapId = "APORNAP-CPLX[c]"
		rnapIndex = moleculeIds.index(rnapId)
		rnapCountsBulk = bulkMolecules.readColumn("counts")[:, rnapIndex]

		bulkMolecules.close()

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		rnapIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		nActive = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnapIndex]

		uniqueMoleculeCounts.close()

		plt.figure(figsize = (8.5, 11))

		plt.plot(time / 60., nActive*100. / ( nActive + rnapCountsBulk))
		#plt.axis([0,60,0,25])
		plt.xlabel("Time (min)")
		plt.ylabel("Percent of RNA Polymerase Molecules that are Active")
		plt.title("Active RNA Polymerase Percentage")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
