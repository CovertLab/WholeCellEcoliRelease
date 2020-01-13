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


		proteinIds = sim_data.process.translation.monomerData["id"]
		rnaIds = sim_data.process.translation.monomerData["rnaId"]

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")
		moleculeIds = bulkMolecules.readAttribute("objectNames")

		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		all_mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')

		rnaIndexes = np.array([all_mRNA_ids.index(moleculeId) for moleculeId in rnaIds], np.int)
		rnaCountsBulk = mRNA_counts[:, rnaIndexes]

		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)
		proteinCountsBulk = bulkMoleculeCounts[:, proteinIndexes]

		relativeMRnaCounts = rnaCountsBulk[-1, :]
		relativeProteinCounts = proteinCountsBulk[-1, :]

		plt.figure(figsize = (8.5, 11))

		plt.plot(relativeMRnaCounts, relativeProteinCounts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("RNA count (at final time step)")
		plt.ylabel("Protein count (at final time step)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
