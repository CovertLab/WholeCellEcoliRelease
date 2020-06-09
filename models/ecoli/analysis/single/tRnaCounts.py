"""
Plot tRNA counts

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/8/2014
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get the names of rnas from the KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		isTRna = sim_data.process.transcription.rnaData["isTRna"]
		rnaIds = sim_data.process.transcription.rnaData["id"][isTRna]
		charged_trna_ids = sim_data.process.transcription.charged_trna_names

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mol_indices = {mol: i for i, mol in enumerate(moleculeIds)}

		uncharged_indices = np.array([mol_indices[moleculeId] for moleculeId in rnaIds], np.int)
		charged_indices = np.array([mol_indices[moleculeId] for moleculeId in charged_trna_ids], np.int)

		bulk_counts = bulkMolecules.readColumn("counts")
		rna_counts = bulk_counts[:, uncharged_indices] + bulk_counts[:, charged_indices]

		plt.figure(figsize = (8.5, 11))

		counts = rna_counts[-1, :]
		expectedCountsArbitrary = sim_data.process.transcription.rnaExpression[sim_data.condition][isTRna]
		expectedCounts = expectedCountsArbitrary/expectedCountsArbitrary.sum() * counts.sum()

		maxLine = 1.1 * max(expectedCounts.max(), counts.max())
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("Expected tRNA count (scaled to total)")
		plt.ylabel("Actual tRNA count (at final time step)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
