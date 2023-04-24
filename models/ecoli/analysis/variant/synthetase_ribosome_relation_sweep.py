"""
Plots relationship between HisRS and ribosome elongation rate when using
different HisRS kcats.
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

		variants_simulated = self.ap.get_variants()
		with open(self.ap.get_variant_kb(variants_simulated[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		variants = [3, 4, 5, 0]
		if not np.all([variant in variants_simulated for variant in variants]):
			print('This analysis script is designed for the variant:'
				f' trna_synthetase_kinetics variant indexes: {variants}.')
			return

		# Get data and plot
		synthetase = 'HISS-CPLX[c]'
		variant_to_data = {}
		variant_to_kcat = {}
		for variant in variants:

			# Get kcat
			with open(self.ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)
			k_cat = (variant_sim_data.relation.synthetase_to_k_cat[synthetase]
				).asNumber(1/units.s)
			variant_to_kcat[variant] = k_cat

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
		color_sweep_low = '#b5cf6b'
		color_sweep_high = '#e7ba52'
		fontsize_legend = 7
		fontsize_label = 9
		alpha_histogram = 0.5
		alpha_scatter = 0.1

		################################################################
		# Scatter: ribosome rate vs synthetase concentration
		variants = [3, 4, 5, 0]
		colors = [color_experiment, color_sweep_low, color_sweep_high,
			color_control]

		ax = axes[1, 0]
		for variant, color in zip(variants, colors):
			x, y = variant_to_data[variant]
			ax.scatter(x, y, marker='.', s=4, linewidths=0, color=color,
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
			mlines.Line2D([], [], color=color_control,
				label=f'k_cat = {round(variant_to_kcat[0])} 1/s (Optimized)',
				**legend_format),
			mlines.Line2D([], [], color=color_sweep_high,
				label=f'k_cat = {round(variant_to_kcat[5])} 1/s',
				**legend_format),
			mlines.Line2D([], [], color=color_sweep_low,
				label=f'k_cat = {round(variant_to_kcat[4])} 1/s',
				**legend_format),
			mlines.Line2D([], [], color=color_experiment,
				label=f'k_cat = {round(variant_to_kcat[3])} 1/s (Measured)',
				**legend_format),
			], loc='lower right', fontsize=7)

		################################################################
		variants = [0, 3, 4, 5]
		props = [
			# Variant 0
			{
			'color': color_control,
			'linewidth': 0,
			'alpha': alpha_histogram,
			},

			# Variant 3
			{
			'edgecolor': color_experiment,
			'linewidth': 1.5,
			'fill': False,
			},

			# Variant 4
			{
			'edgecolor': color_sweep_low,
			'linewidth': 1.5,
			'fill': False,
			},

			# Variant 5
			{
			'edgecolor': color_sweep_high,
			'linewidth': 1.5,
			'fill': False,
			},
		]

		################################################################
		# Histogram: synthetase concentration
		x_min = min([variant_to_data[variant][0].min() for variant in variants])
		x_max = max([variant_to_data[variant][0].max() for variant in variants])
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
		x_min = min([variant_to_data[variant][1].min() for variant in variants])
		x_max = max([variant_to_data[variant][1].max() for variant in variants])
		bins = np.arange(np.floor(x_min), np.ceil(x_max) + 1, 1)

		ax = axes[1, 1]
		ax.set_xlabel('Frequency', fontsize=fontsize_label)
		for variant, prop in zip(variants, props):
			data = variant_to_data[variant][1]
			weights = np.ones_like(data) / len(data)
			ax.hist(data, bins=bins, weights=weights, orientation='horizontal',
				**prop,)
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
		control_min_ribosome_rate = round(
			variant_to_data[0][1].min(), 1)

		for variant in [3, 4, 5, 0]:
			print(f'Variant {variant}')

			kcat = variant_to_kcat[variant]
			print(f'\tkcat:\t{kcat:.1f}')

			# Minimum elongation rate
			y_min = variant_to_data[variant][1].min()
			print(f'\ty min:\t\t{y_min:.1f}')

			if variant != 0:
				mask = variant_to_data[variant][1] < 17
				x_threshold = variant_to_data[variant][0][mask].max()
				print(f'\tx threshold:\t{x_threshold:.2f}')

				# Portion of time steps with ribosome elongation rate <
				# control's min ribosome rate
				mask = variant_to_data[variant][1] < control_min_ribosome_rate
				f = sum(mask) / mask.shape[0] * 100
				print(f'\t<{control_min_ribosome_rate}:\t\t{f:.1f} %')


if __name__ == "__main__":
	Plot().cli()
