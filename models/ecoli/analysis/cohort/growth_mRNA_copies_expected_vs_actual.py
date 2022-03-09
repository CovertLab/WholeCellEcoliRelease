"""
Generates a scatter plot of growth-associated mRNA cistron copy numbers (genes
 that encode for ribosomal proteins and RNAPs) that are expected from the
expression levels calculated in the ParCa vs actual copies in the simulations.
Both values are normalized such that the sum of copy numbers of all plotted mRNA
species sum to one.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader


FIGSIZE = (4, 4)
BOUNDS = [0.5, 2.5]

EXPECTED_COUNT_CONDITION = 'basal'


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		simOutDir = os.path.join(cell_paths[0], 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Calculate expected copies from sim_data
		is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		all_cistron_ids = sim_data.process.transcription.cistron_data['id']
		all_mRNA_cistron_ids = all_cistron_ids[is_mRNA]
		mRNA_cistron_id_to_index = {
			cistron_id: i for i, cistron_id in enumerate(all_mRNA_cistron_ids)
			}
		rprotein_cistron_ids = all_cistron_ids[
			sim_data.process.transcription.cistron_data['is_ribosomal_protein']]
		rnap_cistron_ids = all_cistron_ids[
			sim_data.process.transcription.cistron_data['is_RNAP']]
		rprotein_indexes = np.array(
			[mRNA_cistron_id_to_index[x] for x in rprotein_cistron_ids])
		rnap_indexes = np.array(
			[mRNA_cistron_id_to_index[x] for x in rnap_cistron_ids])

		expected_mRNA_counts = sim_data.process.transcription.fit_cistron_expression[EXPECTED_COUNT_CONDITION][is_mRNA]
		expected_rprotein_counts = expected_mRNA_counts[rprotein_indexes]
		expected_rnap_counts = expected_mRNA_counts[rnap_indexes]

		assert np.all(
			mRNA_cistron_ids == sim_data.process.transcription.cistron_data['id'][is_mRNA])

		# Read actual counts
		all_actual_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts')
		all_actual_rprotein_counts = all_actual_counts[:, rprotein_indexes]
		all_actual_rnap_counts = all_actual_counts[:, rnap_indexes]

		# Get average count across all timesteps across all sims
		actual_rprotein_counts = all_actual_rprotein_counts.mean(axis=0)
		actual_rnap_counts = all_actual_rnap_counts.mean(axis=0)

		# Normalize expected copies to total sum of actual counts
		expected_sum = expected_rprotein_counts.sum() + expected_rnap_counts.sum()
		actual_sum = actual_rprotein_counts.sum() + actual_rnap_counts.sum()
		expected_rprotein_counts *= actual_sum/expected_sum
		expected_rnap_counts *= actual_sum/expected_sum

		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)

		ax.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
		ax.scatter(
			np.log10(expected_rprotein_counts + 1),
			np.log10(actual_rprotein_counts + 1),
			s=3, c='r',
			label='ribosomal proteins')
		ax.scatter(
			np.log10(expected_rnap_counts + 1),
			np.log10(actual_rnap_counts + 1),
			s=3, c='b',
			label='RNAP subunits')

		ax.set_title('Expected vs actual RNA copies')
		ax.set_xlabel('$\log_{10}$(Expected copies + 1)')
		ax.set_ylabel('$\log_{10}$(Actual copies + 1)')
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_position(("outward", 20))
		ax.spines["left"].set_position(("outward", 20))
		ax.set_xlim(BOUNDS)
		ax.set_ylim(BOUNDS)
		ax.legend(loc=2)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
