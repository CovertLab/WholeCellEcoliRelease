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

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader



FIGSIZE = (6, 6)
BOUNDS = [1e-9, 1]
NUMERICAL_ZERO = 1e-10

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
		assert np.all(
			mRNA_cistron_ids == sim_data.process.transcription.cistron_data['id'][is_mRNA])
		expected_counts = sim_data.process.transcription.fit_cistron_expression[EXPECTED_COUNT_CONDITION][is_mRNA]

		# Normalize counts
		expected_counts /= expected_counts.sum()

		# Read actual counts
		all_actual_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts', fun=lambda x: x.mean(axis=0))

		# Get average count across all sims
		actual_counts = all_actual_counts.mean(axis=0)

		# Normalize counts
		actual_counts /= actual_counts.sum()

		plt.figure(figsize=FIGSIZE)

		plt.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
		plt.scatter(
			expected_counts + NUMERICAL_ZERO,
			actual_counts + NUMERICAL_ZERO,
			alpha=0.5,
			s=1)

		plt.title('Expected vs actual RNA copies')
		plt.xlabel('Expected normalized copies')
		plt.ylabel('Actual normalized copies')
		plt.xlim(BOUNDS)
		plt.ylim(BOUNDS)
		plt.xscale('log')
		plt.yscale('log')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
