"""
Plot mRNA counts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


# TODO: account for complexation

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get the names of rnas from the KB

		sim_data = cPickle.load(open(simDataFile, "rb"))

		isMRna = sim_data.process.transcription.rnaData["isMRna"]

		rnaIds = sim_data.process.transcription.rnaData["id"][isMRna]

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		moleculeIds = bulkMolecules.readAttribute("objectNames")

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]

		# avgCounts = rnaCountsBulk.mean(0)

		# relativeCounts = avgCounts / avgCounts.sum()

		# relativeCounts = rnaCountsBulk[-1, :] / rnaCountsBulk[-1, :].sum()

		bulkMolecules.close()

		plt.figure(figsize = (8.5, 11))

		counts = rnaCountsBulk[-1, :]

		expectedCountsArbitrary = sim_data.process.transcription.rnaExpression[sim_data.condition][isMRna]

		expectedCounts = expectedCountsArbitrary/expectedCountsArbitrary.sum() * counts.sum()

		maxLine = 1.1 * max(expectedCounts.max(), counts.max())
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("Expected RNA count (scaled to total)")
		plt.ylabel("Actual RNA count (at final time step)")

		# plt.show()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
