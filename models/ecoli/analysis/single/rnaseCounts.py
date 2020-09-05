"""
Plot RNAse counts

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/14/2015
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, "rb"))

		endoRnaseIds = sim_data.process.rna_decay.endoRNase_ids
		exoRnaseIds = sim_data.molecule_groups.exoRNases
		RNase_IDS = np.concatenate((endoRnaseIds, exoRnaseIds))

		endoRnase_RnaIDs = sim_data.molecule_groups.endoRNase_rnas
		exoRnase_RnaIDs = sim_data.molecule_groups.exoRNase_rnas
		RNase_RnaIDS = np.concatenate((endoRnase_RnaIDs, exoRnase_RnaIDs))

		# Load count data for mRNAs
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		all_mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')

		# Get counts for RNase proteins and mRNAs
		(RNase_counts,) = read_bulk_molecule_counts(simOutDir, (RNase_IDS,))
		rnaIndexes = np.array([all_mRNA_ids.index(rna) for rna in RNase_RnaIDS], np.int)
		RNase_RNA_counts = mRNA_counts[:, rnaIndexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		n_subplots = 2*len(RNase_IDS)

		plt.figure(figsize = (8.5, 11))
		plt.rc('xtick', labelsize=7)
		plt.rc('ytick', labelsize=5)

		for subplotIdx in range(0, n_subplots):
			ax = plt.subplot(18, 2, 1 + subplotIdx)

			if not subplotIdx % 2:
				plt.plot(time / 60., RNase_counts[:, subplotIdx // 2])
			else:
				plt.plot(time / 60., RNase_RNA_counts[:, subplotIdx // 2])

			if not subplotIdx >= n_subplots - 2:
				frame = plt.gca()
				for xlabel_i in frame.axes.get_xticklines():
					xlabel_i.set_visible(True)
				for xlabel_i in frame.axes.get_xticklabels():
					xlabel_i.set_visible(False)

			if subplotIdx >= n_subplots - 2:
				plt.xlabel("Time (min)", fontsize = 7)

			if not subplotIdx % 2:
				plt.ylabel("Protein counts", fontsize = 5)
				plt.title(RNase_IDS[subplotIdx // 2], fontsize=7)
			else:
				plt.ylabel("RNA counts", fontsize = 5)
				plt.title(RNase_RnaIDS[subplotIdx // 2], fontsize=7)

			max_yticks = 4
			yloc = plt.MaxNLocator(max_yticks)
			ax.yaxis.set_major_locator(yloc)

		plt.subplots_adjust(hspace = 0.75, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
