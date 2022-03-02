"""
Generates a scatter plot of mRNA cistron copy numbers that are expected from the
expression levels calculated in the ParCa vs actual copies in the simulations.
Both values are normalized such that the sum of copy numbers of all mRNA species
sum to one.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import scipy.stats as st

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader


FIGSIZE = (6, 6)
BOUNDS = [-0.5, 3]
P_VALUE_THRESHOLD = 1e-3

EXPECTED_COUNT_CONDITION = 'basal'


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		simOutDir = os.path.join(cell_paths[0], 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Get mask for mRNA genes that are
		# i) Does not encode for ribosomal proteins or RNAPs
		# ii) not affected by manual overexpression/underexpression
		# to isolate the effects of operons on expression levels.
		is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		assert np.all(
			mRNA_cistron_ids == sim_data.process.transcription.cistron_data['id'][is_mRNA])

		mRNA_is_rnap_or_rprotein = np.logical_or(
				sim_data.process.transcription.cistron_data['is_RNAP'],
				sim_data.process.transcription.cistron_data['is_ribosomal_protein'])[is_mRNA]

		is_adjusted = np.zeros_like(is_mRNA, dtype=bool)
		all_rna_ids = sim_data.process.transcription.rna_data['id']
		for adjusted_cistron_id in sim_data.adjustments.rna_expression_adjustments.keys():
			# Include cistrons whose expression is adjusted because they belong
			# to the same TU as the cistron that is bumped up
			adjusted_rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(
				adjusted_cistron_id)
			adjusted_cistron_indexes = []
			for adjusted_rna_index in adjusted_rna_indexes:
				adjusted_cistron_indexes.extend(
					sim_data.process.transcription.rna_id_to_cistron_indexes(
						all_rna_ids[adjusted_rna_index]))

			is_adjusted[adjusted_cistron_indexes] = True

		mRNA_is_adjusted = is_adjusted[is_mRNA]
		mRNA_mask = np.logical_and(~mRNA_is_rnap_or_rprotein, ~mRNA_is_adjusted)

		# Get list of cistron IDs that are plotted
		plotted_mRNA_cistrons = []
		for i in np.where(mRNA_mask)[0]:
			plotted_mRNA_cistrons.append(mRNA_cistron_ids[i])

		# Get expected counts
		expected_counts = sim_data.process.transcription.cistron_expression[EXPECTED_COUNT_CONDITION][
			is_mRNA][mRNA_mask]

		# Read actual counts
		all_actual_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts')

		# Get average count across all timesteps over all sims
		actual_counts = all_actual_counts[:, mRNA_mask].mean(axis=0)
		n_timesteps = all_actual_counts.shape[0]

		# Normalize expected counts with sum of actual counts
		expected_counts = expected_counts/expected_counts.sum() * actual_counts.sum()

		# Get statistical significance boundaries assuming a Poissonian
		# distribution
		z_score_threshold = st.norm.ppf(1 - P_VALUE_THRESHOLD/2)
		ub = expected_counts + z_score_threshold*np.sqrt(expected_counts/n_timesteps)
		lb = expected_counts - z_score_threshold*np.sqrt(expected_counts/n_timesteps)

		# Get mask for outliers
		outlier_mask = np.logical_or(actual_counts > ub, actual_counts < lb)

		plt.figure(figsize=FIGSIZE)

		plt.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
		plt.scatter(
			np.log10(expected_counts[~outlier_mask] + 1),
			np.log10(actual_counts[~outlier_mask] + 1),
			c='#cccccc', s=1, label=f'p ≥ {P_VALUE_THRESHOLD:g}')
		# Highlight outliers in blue
		plt.scatter(
			np.log10(expected_counts[outlier_mask] + 1),
			np.log10(actual_counts[outlier_mask] + 1),
			c='b', s=1, label=f'p < {P_VALUE_THRESHOLD:g}')

		plt.title('Expected vs actual RNA copies')
		plt.xlabel('Expected normalized copies')
		plt.ylabel('Actual normalized copies')
		plt.xlim(BOUNDS)
		plt.ylim(BOUNDS)
		plt.legend()

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
