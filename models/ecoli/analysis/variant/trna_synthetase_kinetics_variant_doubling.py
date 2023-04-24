"""
Compares doubling times of trna synthetase kinetics variants.
"""

import pickle
import os
import re

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

print_stats = True

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# Check for the required variants
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		if not (0 in variants and 1 in variants and 2 in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0, 1, and 2.')
			return

		# Get data
		variant_control = 0
		variant_previous = 1
		variant_measured = 2

		variants = [variant_control, variant_previous, variant_measured]
		variant_to_times = {}

		for variant in variants:

			# Get cell paths
			cell_paths = ap.get_cells(variant=[variant])

			# Doubling times
			variant_to_times[variant] = []
			for sim_dir in cell_paths:
				sim_out_dir = os.path.join(sim_dir, 'simOut')
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				time = main_reader.readColumn('time') / 60
				variant_to_times[variant].append(time[-1] - time[0])

		guide = sim_data.condition_to_doubling_time[sim_data.condition]
		guide = guide.asNumber(units.min)

		# Calculate z-score: (mean - population_mean) / std
		population_mean = 44 # measurement

		print('z scores')
		for variant in variants:
			sample = variant_to_times[variant]
			mean = np.mean(sample)
			std = np.std(sample)

			z = (mean - population_mean) / std
			print(f'\tVariant {variant}:\t{z:.2f}')

		# Plot
		fig, axes = plt.subplots(1, 2, figsize=(5, 2.5))
		color_control = '#bdbdbd'
		color_experiment = '#6b6ecf'
		alpha_control = 0.7
		alpha_experiment = 0.6
		fontsize_heading = 11
		fontsize_label = 9
		fontsize_ticks = 7

		hist_params = {
			'edgecolor': 'w',
			'linewidth': 0.2,
		}

		variant_to_weights = {}
		for variant in variants:
			weights = np.ones_like(variant_to_times[variant])
			weights /= weights.shape[0]
			variant_to_weights[variant] = weights

		################################################################
		# Figure 2
		ax = axes[0]

		# Previous
		bin_size = 3
		bins = np.arange(
			np.floor(min(variant_to_times[variant_previous]) / 10) * 10,
			np.ceil(max(variant_to_times[variant_previous]) / 10) * 10 + bin_size,
			bin_size)

		ax.hist(variant_to_times[variant_previous],
			weights=variant_to_weights[variant_previous],
			bins=bins,
			color=color_control,
			label='Previous Model',
			**hist_params,
			)

		# Measured
		bin_size = 10
		bins = np.arange(
			np.floor(min(variant_to_times[variant_measured]) / 10) * 10,
			np.ceil(max(variant_to_times[variant_measured]) / 10) * 10 + bin_size,
			bin_size)
		ax.hist(variant_to_times[variant_measured],
			weights=variant_to_weights[variant_measured],
			bins=bins,
			color=color_experiment,
			label='Measured tRNA Synthetase kcats',
			alpha=alpha_experiment,
			**hist_params,
			)
		ax.axvline(guide, color='k', linewidth=1)
		ax.legend(fontsize=4)
		ax.set_xticks([50, 100, 150, 200])
		ax.set_yticks([0, 0.1, 0.2, 0.3])
		ax.set_ylim([0, 0.32])

		################################################################
		# Figure 3
		ax = axes[1]

		# Previous
		bin_size = 2
		bins = np.arange(
			np.floor(min(variant_to_times[variant_previous]) / 10) * 10,
			np.ceil(max(variant_to_times[variant_previous]) / 10) * 10 + bin_size,
			bin_size)

		ax.hist(variant_to_times[variant_previous],
			weights=variant_to_weights[variant_previous],
			bins=bins,
			color=color_control,
			alpha=0.6, 
			label='Previous Model',
			**hist_params,
			)

		# Current
		bin_size = 2
		bins = np.arange(
			np.floor(min(variant_to_times[variant_control]) / 10) * 10,
			np.ceil(max(variant_to_times[variant_control]) / 10) * 10 + bin_size,
			bin_size)

		ax.hist(variant_to_times[variant_control],
			weights=variant_to_weights[variant_control],
			bins=bins,
			label='Current Model',
			# alpha=alpha_experiment,
			# **hist_params,
			facecolor='none',
			edgecolor=color_experiment,
			linewidth=0.7,
			)
		ax.axvline(guide, color='k', linewidth=1)
		ax.legend(fontsize=4)
		ax.set_yticks([0, 0.1, 0.2])

		################################################################

		for ax in axes:
			ax.set_title('Doubling Time', fontsize=fontsize_heading)
			ax.set_xlabel('Time (min)', fontsize=fontsize_label)
			ax.set_ylabel('Fraction of Simulated Cells', fontsize=fontsize_label)

			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		if print_stats:
			previous = np.mean(variant_to_times[variant_previous])
			std = np.std(variant_to_times[variant_previous])
			print(f'Previous: {previous:.2f}')
			print(f'\tstd: {std:.2f}')

			measured = np.array(variant_to_times[variant_measured])
			f_at_180 = sum(measured == 180) / measured.shape[0]
			print(f'Measured kcats: {measured.mean():.2f} -- {measured.mean() / previous:.2f} fold')
			print(f'\tPercent at 180: {f_at_180 * 100:.3f}')

			std = np.std(variant_to_times[variant_control])
			print(f'Current: {np.mean(variant_to_times[variant_control]):.2f}')
			print(f'\tstd: {std:.2f}')


if __name__ == "__main__":
	Plot().cli()
