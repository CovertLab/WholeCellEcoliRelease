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

		rnaIds = sim_data.process.transcription.rnaData["id"][sim_data.relation.rnaIndexToMonomerMapping]

		proteinIds = sim_data.process.translation.monomerData["id"]

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")

		moleculeIds = bulkMolecules.readAttribute("objectNames")

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMoleculeCounts[:, rnaIndexes]

		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

		proteinCountsBulk = bulkMoleculeCounts[:, proteinIndexes]

		bulkMolecules.close()

		relativeMRnaCounts = rnaCountsBulk[-1, :] #/ rnaCountsBulk[-1, :].sum()
		relativeProteinCounts = proteinCountsBulk[-1, :] #/ proteinCountsBulk[-1, :].sum()

		plt.figure(figsize = (8.5, 11))

		plt.plot(relativeMRnaCounts, relativeProteinCounts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("RNA count (at final time step)")
		plt.ylabel("Protein count (at final time step)")

		# plt.show()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
