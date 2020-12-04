"""
Compare fold changes from different sources.
"""

from six.moves import cPickle
import os

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = cPickle.load(f)

		# Load original fold change data
		original_fold_changes = {}
		consistent = {}
		for fc in raw_data.fold_changes:
			pair = (fc['TF'], fc['Target'])
			original_fold_changes[pair] = fc['log2 FC mean']
			consistent[pair] = fc['Regulation_direct'] <= 2

		# Load NCA fold change data
		nca_fold_changes = {}
		for fc in raw_data.fold_changes_nca:
			nca_fold_changes[(fc['TF'], fc['Target'])] = fc['log2 FC mean']

		# Compare regulation pairs in both datasets
		original = []
		nca = []
		for pair in original_fold_changes:
			if pair in nca_fold_changes and consistent[pair]:
				original.append(original_fold_changes[pair])
				nca.append(nca_fold_changes[pair])
		original = np.array(original)
		nca = np.array(nca)
		pearson = stats.pearsonr(original, nca)

		# Get data unique to each dataset
		original_unique_pairs = set(original_fold_changes) - set(nca_fold_changes)
		nca_unique_pairs = set(nca_fold_changes) - set(original_fold_changes)
		original_unique = np.array([original_fold_changes[pair] for pair in original_unique_pairs])
		nca_unique = np.array([nca_fold_changes[pair] for pair in nca_unique_pairs])

		# Get range for y=x line
		min_val = min(original.min(), nca.min())
		max_val = max(original.max(), nca.max())
		unique_min = min(original_unique.min(), nca_unique.min())
		unique_max = max(original_unique.max(), nca_unique.max())

		plt.figure(figsize=(5, 15))
		gs = gridspec.GridSpec(3, 1)
		scatter_ax = plt.subplot(gs[0, 0])
		shared_hist_ax = plt.subplot(gs[1, 0])
		unique_hist_ax = plt.subplot(gs[2, 0])

		# Scatter plot
		## Plot comparison
		scatter_ax.plot(original, nca, 'x', label='log2 fold changes')
		xlim = scatter_ax.get_xlim()
		ylim = scatter_ax.get_ylim()

		## Plot y=x and revert to previous limits
		scatter_ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='y=x')
		scatter_ax.set_xlim(xlim)
		scatter_ax.set_ylim(ylim)

		## Format axes
		scatter_ax.set_xlabel('Original fold change')
		scatter_ax.set_ylabel('NCA predicted fold change')
		scatter_ax.spines['top'].set_visible(False)
		scatter_ax.spines['right'].set_visible(False)
		scatter_ax.legend(fontsize=8, frameon=False)

		## Display stats
		scatter_ax.set_title(f'Pearson r={pearson[0]:.3f} (p={pearson[1]:.0e}, n={len(original)})', fontsize=10)

		# Histogram for shared fold change pairs
		## Plot data
		shared_hist_ax.hist(original, range=(min_val, max_val), bins=20, alpha=0.5, label='Shared original fold changes')
		shared_hist_ax.hist(nca, range=(min_val, max_val), bins=20, alpha=0.5, label='Shared NCA fold changes')

		## Format axes
		shared_hist_ax.set_xlabel('Fold change (shared regulation in datasets)')
		shared_hist_ax.set_ylabel('# of regulatory pairs')
		shared_hist_ax.spines['top'].set_visible(False)
		shared_hist_ax.spines['right'].set_visible(False)
		shared_hist_ax.legend(fontsize=8, frameon=False)

		# Histogram for shared fold change pairs
		## Plot data
		unique_hist_ax.hist(original_unique, range=(unique_min, unique_max),
			bins=20, alpha=0.5, label=f'Unique original fold changes (N={len(original_unique)})')
		unique_hist_ax.hist(nca_unique, range=(unique_min, unique_max),
			bins=20, alpha=0.5, label=f'Unique NCA fold changes (N={len(nca_unique)})')

		## Format axes
		unique_hist_ax.set_xlabel('Fold change (unique regulation in datasets)')
		unique_hist_ax.set_ylabel('# of regulatory pairs')
		unique_hist_ax.spines['top'].set_visible(False)
		unique_hist_ax.spines['right'].set_visible(False)
		unique_hist_ax.legend(fontsize=8, frameon=False)

		# Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
