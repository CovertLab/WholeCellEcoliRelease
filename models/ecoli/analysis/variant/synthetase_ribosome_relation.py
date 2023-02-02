"""
Plots relationship between HisRS and ribosome elongation rate.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.lines as mlines
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

print_stats = True

def get_fit(x, k, theta):
	numerator = np.power(x, k-1) * np.exp(-x / theta)
	denominator = gamma_function(k) * (theta**k)
	return numerator / denominator

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		variants = self.ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		variant_control = 0
		variant_experiment = 3
		# Analyze variants 0 (control) and 3 (HisRS drop)
		if not (variant_control in variants and variant_experiment in variants):
			print('This analysis script is designed for the variant:'
				f' trna_synthetase_kinetics variant indexes {variant_control}'
				f' and {variant_experiment}.')
			return

		# Get data and plot
		synthetase = 'HISS-CPLX[c]'
		variant_to_data = {}
		for variant in [variant_control, variant_experiment]:
			cell_paths = self.ap.get_cells(variant=[variant])
			ribosome_rate = read_stacked_columns(
				cell_paths, 'RibosomeData', 'effectiveElongationRate',
				remove_first=True).reshape(-1)
			(n_synthetase,) = read_stacked_bulk_molecules(
				cell_paths, ([synthetase],), remove_first=True)

			cell_volume = 1e-15 * units.L * read_stacked_columns(
				cell_paths, 'Mass', 'cellVolume', remove_first=True).reshape(-1)
			c_synthetase = (1
				/ sim_data.constants.n_avogadro
				/ cell_volume
				* n_synthetase
				).asNumber(units.umol / units.L)
			variant_to_data[variant] = (c_synthetase, ribosome_rate)

		# Plot
		fig, axes = plt.subplots(2, 2, figsize=(3.5, 3),
			gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [3, 1]})

		color_control = '#bdbdbd'
		color_experiment = '#6b6ecf'
		fontsize_legend = 7
		fontsize_label = 9
		alpha_histogram = 0.5
		alpha_scatter = 0.1

		variants = [0, 3]

		control_histogram_props = {
			'color': color_control,
			'linewidth': 0,
			'alpha': alpha_histogram,
		}
		experiment_histogram_props = {
			'edgecolor': color_experiment,
			'linewidth': 1.5,
			'fill': False,
		}

		props = [control_histogram_props, experiment_histogram_props]

		################################################################
		# Scatter: ribosome rate vs synthetase concentration
		ax = axes[1, 0]
		x, y = variant_to_data[3]
		ax.scatter(x, y, marker='.', s=4, linewidths=0, color=color_experiment,
			alpha=alpha_scatter, rasterized=True)

		x, y = variant_to_data[0]
		ax.scatter(x, y, marker='.', s=4, linewidths=0, color=color_control,
			alpha=alpha_scatter, rasterized=True)

		ax.set_xlabel('HisRS Concentration (uM)', fontsize=fontsize_label)
		ax.set_ylabel('Ribosome Elongation Rate\n(amino acids/s/ribosome)',
			fontsize=fontsize_label)

		# Legend
		legend_format = {
			'marker': '.',
			'linestyle': 'None',
			'markersize': 4,
			}
		ax.legend(handles=[
			mlines.Line2D([], [], color=color_control, label='Optimized',
				**legend_format),
			mlines.Line2D([], [], color=color_experiment, label='Measured',
				**legend_format)
			], loc='lower right', fontsize=7)

		################################################################
		# Histogram: synthetase concentration
		x_min = min(
			variant_to_data[variant_control][0].min(),
			variant_to_data[variant_experiment][0].min())
		x_max = max(
			variant_to_data[variant_control][0].max(),
			variant_to_data[variant_experiment][0].max())
		bins = np.arange(np.floor(x_min), np.ceil(x_max) + 0.025, 0.025)

		ax = axes[0, 0]
		ax.set_ylabel('Frequency', fontsize=fontsize_label)
		for variant, prop in zip(variants, props):
			data = variant_to_data[variant][0]
			weights = np.ones_like(data) / len(data)
			ax.hist(data, bins=bins, weights=weights, **prop)
		ax.set_xlim(axes[1, 0].get_xlim())
		ax.set_yticks([0, 0.1, 0.2])
		ax.set_ylim([0, 0.2])

		################################################################
		# Blank
		ax = axes[0, 1]
		ax.remove()

		################################################################
		# Histogram: ribosome rate
		x_min = min(
			variant_to_data[variant_control][1].min(),
			variant_to_data[variant_experiment][1].min())
		x_max = max(
			variant_to_data[variant_control][1].max(),
			variant_to_data[variant_experiment][1].max())
		bins = np.arange(np.floor(x_min), np.ceil(x_max) + 1, 1)

		ax = axes[1, 1]
		ax.set_xlabel('Frequency', fontsize=fontsize_label)
		for variant, prop in zip(variants, props):
			data = variant_to_data[variant][1]
			weights = np.ones_like(data) / len(data)
			ax.hist(data, bins=bins, weights=weights, orientation='horizontal',
				**prop)
		ax.set_ylim(axes[1, 0].get_ylim())

		################################################################
		# Format
		for ax in axes.reshape(-1):
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(axis='both', which='major', labelsize=fontsize_legend)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Summary stats
		if print_stats:
			control_min_ribosome_rate = round(
				variant_to_data[variant_control][1].min(), 1)

			mask = (variant_to_data[variant_experiment][1]
				< control_min_ribosome_rate)
			threshold = variant_to_data[variant_experiment][0][mask].max()
			print(f'HisRS concentration threshold: {threshold:.2f}')
			print(f'Portion of time steps below {control_min_ribosome_rate}'
				f' aa/s: {sum(mask) / len(mask) * 100:.1f}%')

			x_min = variant_to_data[variant_experiment][0][mask].min()
			x_med = variant_to_data[variant_experiment][0][mask].max()
			x_max = variant_to_data[variant_experiment][0].max()
			print(f'HisRS concentration:\n\t{x_min:.1f} to {x_med:.1f}'
				f' ({(x_med - x_min) / (x_max - x_min) * 100:.1f}% of range)')

			y_min = variant_to_data[variant_experiment][1][mask].min()
			y_med = variant_to_data[variant_experiment][1][mask].max()
			y_max = variant_to_data[variant_experiment][1].max()
			print(f'Ribosome elongation rate:\n\t{y_min:.2f} to {y_med:.2f}'
				f' ({(y_med - y_min) / (y_max - y_min) * 100:.1f}% of range)')


if __name__ == "__main__":
	Plot().cli()
