"""
Plot mRNA counts
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


		protein_ids = sim_data.process.translation.monomer_data["id"]
		cistron_ids = sim_data.process.translation.monomer_data['cistron_id']

		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_counts = mRNA_counts_reader.readColumn('mRNA_cistron_counts')
		all_mRNA_cistron_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}

		rna_cistron_indexes = np.array([all_mRNA_cistron_idx[moleculeId] for moleculeId in cistron_ids], int)
		rna_cistron_counts_bulk = mRNA_cistron_counts[:, rna_cistron_indexes]

		(protein_counts_bulk,) = read_bulk_molecule_counts(simOutDir, (protein_ids,))

		relativeMRnaCounts = rna_cistron_counts_bulk[-1, :]
		relativeProteinCounts = protein_counts_bulk[-1, :]

		plt.figure(figsize = (8.5, 11))

		plt.plot(relativeMRnaCounts, relativeProteinCounts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("RNA count (at final time step)")
		plt.ylabel("Protein count (at final time step)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
