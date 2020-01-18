'''
Patterns of co-transcriptional translation in the E. coli WCM

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/17/20
'''

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

PLOT_TOP_N_GENES = 30


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')
		all_mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		full_mRNA_counts = mRNA_counts_reader.readColumn('full_mRNA_counts')
		partial_mRNA_counts = mRNA_counts_reader.readColumn('partial_mRNA_counts')

		fig = plt.figure()
		fig.set_size_inches(8, 11.5)

		gs = gridspec.GridSpec(2, 1)

		# Plot counts of full/partial transcripts over time
		ax = plt.subplot(gs[0, 0])
		ax.plot(time, all_mRNA_counts.sum(axis=1), label='All mRNAs')
		ax.plot(time, full_mRNA_counts.sum(axis=1), label='Fully transcribed mRNAs')
		ax.plot(time, partial_mRNA_counts.sum(axis=1), label='Partially transcribed mRNAs')
		ax.legend()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylabel('Total Counts')

		# Get top N genes with the highest percentage of active mRNAs that are
		# partially transcribed (minimum 10 average total transcripts)
		all_mRNA_counts_mean = all_mRNA_counts.mean(axis=0)
		partial_mRNA_proportions = np.nan_to_num(
			partial_mRNA_counts.mean(axis=0)/all_mRNA_counts_mean)
		partial_mRNA_proportions[all_mRNA_counts_mean < 10] = 0
		rank = np.argsort(partial_mRNA_proportions)[::-1][:PLOT_TOP_N_GENES]
		ranked_genes = [sim_data.fathom.common_names[mRNA_ids[i]][0][:-5] for i in rank]

		ax = plt.subplot(gs[1, 0])
		ax.bar(range(PLOT_TOP_N_GENES), partial_mRNA_proportions[rank])
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(ranked_genes, rotation=90)
		ax.set_title('Genes with highest percentage of partially transcribed active transcripts\n($\geq 10$ total transcripts)')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylim([0, 1])
		ax.set_ylabel('Percentage of partially transcribed transcripts')

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
