"""
Plots doubling times as a beeswarm plot.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


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

		# Get doubling times
		variants = [variant_control, variant_experiment]
		data = {variant: [] for variant in variants}
		seeds = []
		sim_dirs = []
		for variant in variants:
			cell_paths = self.ap.get_cells(variant=[variant])

			for sim_dir in cell_paths:
				main_reader = TableReader(
					os.path.join(sim_dir, 'simOut', 'Main'))
				time = main_reader.readColumn('time') / 60
				data[variant].append(time[-1] - time[0])

				if variant == variant_experiment:
					seeds.append(int(sim_dir.split('/')[-3]))
					sim_dirs.append(sim_dir)

		print(f'Top outlier: {max(data[variant_experiment]):.1f}')

		# Beeswarm
		def beeswarm(data, center):
			# Round to nearest full minute
			data = np.round(data).astype(np.int64)

			# Count
			value_to_count = {}
			for i, value in enumerate(np.bincount(data)):
				if value != 0:
					value_to_count[i] = value

			# Assign coordinates
			xs = np.array([])
			ys = np.array([])
			space = 0.05

			for y, count in value_to_count.items():

				# Even
				if count % 2 == 0:
					left = center - space * ((count / 2) - 1) - (space / 2)
					right = left + (space * (count - 1))

					xs = np.hstack((xs, np.linspace(left, right, count)))
					ys = np.hstack((ys, y * np.ones(count)))

				# Odd
				else:
					left = center - (space * ((count - 1) / 2))
					right = left + (space * (count - 1))

					xs = np.hstack((xs, np.linspace(left, right, count)))
					ys = np.hstack((ys, y * np.ones(count)))

			return xs, ys

		def plot_distributions(ax, data, center, color):
			mean = np.mean(data)
			std = np.std(data)
			print(f'\tmean = {mean:.1f}\tstd = {std:.1f}')
			size = 0.6
			z = 3
			ax.plot([center - size, center + size], [mean, mean], color=color,
				linewidth=0.5)
			ax.plot([center - size / 2, center + size / 2],
				[mean + (z * std), mean + (z * std)], color=color,
				linewidth=0.5)
			ax.plot([center - size / 2, center + size / 2],
				[mean - (z * std), mean - (z * std)], color=color,
				linewidth=0.5)
			return

		# Plot
		fig, ax = plt.subplots(1, 1, figsize=(4.5, 2.5))
		color_control = '#c7c7c7'
		color_experiment = '#6b6ecf'
		color_control_accent = '#636363'
		color_experiment_accent = '#393b79'

		xs, ys = beeswarm(data[variant_control], 0)
		ax.scatter(xs, ys, color=color_control, s=5, linewidths=0)
		print('Control')
		plot_distributions(ax, data[variant_control], 0, color_control_accent)

		xs, ys = beeswarm(data[variant_experiment], 1.5)
		ax.scatter(xs, ys, color=color_experiment, s=5, linewidths=0)
		print('Experiment')
		plot_distributions(ax, data[variant_experiment], 1.5,
			color_experiment_accent)

		ax.set_xticks([0, 1.5])
		ax.set_xlim([-0.75, 2.25])
		ax.set_yticks([30, 40, 50, 60, 70])
		ax.set_title('Doubling Time', fontsize=11)
		ax.set_xticklabels(['Optimized', 'Measured'], fontsize=7)
		ax.set_xlabel('HisRS kcat', fontsize=9)
		ax.set_ylabel('Time (min)', fontsize=9)
		ax.tick_params(axis='both', which='major', labelsize=7)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
