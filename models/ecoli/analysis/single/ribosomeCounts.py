"""
Plot ribosome counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
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

		# bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		#
		# moleculeIds = bulkMolecules.readAttribute("objectNames")
		#
		# RIBOSOME_RNA_IDS = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]
		# ribosomeRnaIndexes = np.array([moleculeIds.index(rRnaId) for rRnaId in RIBOSOME_RNA_IDS], np.int)
		# ribosomeRnaCountsBulk = bulkMolecules.readColumn("counts")[:, ribosomeRnaIndexes]

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		nActive = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		plt.figure(figsize = (8.5, 11))

		plt.plot(time / 60, nActive)
		plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("Counts")
		plt.title("Active Ribosomes Final:Initial = %0.2f" % (nActive[-1] / float(nActive[0])))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
