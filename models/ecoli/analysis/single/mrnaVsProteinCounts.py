"""
Plot mRNA counts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


# TODO: account for complexation


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get the names of rnas from the KB

		sim_data = cPickle.load(open(simDataFile, "rb"))


		proteinIds = sim_data.process.translation.monomer_data["id"]
		rnaIds = sim_data.process.translation.monomer_data['rna_id']

		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		all_mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}

		rnaIndexes = np.array([all_mRNA_idx[moleculeId] for moleculeId in rnaIds], np.int)
		rnaCountsBulk = mRNA_counts[:, rnaIndexes]

		(proteinCountsBulk,) = read_bulk_molecule_counts(simOutDir, (proteinIds,))

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
