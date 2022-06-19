"""
Plot RNA polymerase counts and counts of mRNA precursors
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import range


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		rnapId = ["APORNAP-CPLX[c]"]
		(rnapCountsBulk,) = read_bulk_molecule_counts(simOutDir, (rnapId,))

		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_counts = mRNA_counts_reader.readColumn('mRNA_cistron_counts')
		mRNA_cistron_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}

		RNAP_cistron_IDS = [
			"EG10893_RNA", "EG10894_RNA",
			"EG10895_RNA", "EG10896_RNA"]

		rnap_cistron_indexes = np.array([mRNA_cistron_idx[rnapRnaId] for rnapRnaId in RNAP_cistron_IDS], int)
		rnap_cistron_counts = mRNA_cistron_counts[:, rnap_cistron_indexes]

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		rnapIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_RNAP')

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

		for subplotIdx in range(2, 6):
			rnapRnaCountsIdx = subplotIdx - 2

			plt.subplot(5, 1, subplotIdx)

			plt.plot(time / 60., rnap_cistron_counts[:, rnapRnaCountsIdx])
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA counts")
			plt.title(RNAP_cistron_IDS[rnapRnaCountsIdx])

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
