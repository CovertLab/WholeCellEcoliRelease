"""
Plots fold change of mass across cell components, and compares with a control.
"""

import pickle
import os

from matplotlib import pyplot as plt
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

		variant_control = 0
		variant_experiment = 6

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		if not (variant_control in variants and variant_experiment in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 6.')
			return

		# Initialize plot
		fig, axes = plt.subplots(
			3, 1,
			figsize=(3, 2.5),
			sharex='col',
			)

		error_kw = {'lw': 0.5, 'capsize': 2, 'capthick': 0.5}

		component_to_label = {
			'dnaMass': 'DNA',
			'rnaMass': 'RNA',
			'proteinMass': 'Protein',
			}
		components = list(component_to_label.keys())

		def get_bar_parameters(data, per_generation=True):

			if per_generation:
				median = np.nanmedian(data, axis=0)
				q1 = np.nanquantile(data, 0.25, axis=0)
				q3 = np.nanquantile(data, 0.75, axis=0)
				yerr = np.array([median - q1, q3 - median])
			else:
				median = np.nanmedian(data)
				q1 = np.nanquantile(data, 0.25)
				q3 = np.nanquantile(data, 0.75)
				yerr = np.array([median - q1, q3 - median]).reshape((2, 1))

			return median, yerr

		################################################################
		# Variant 0, control
		color_control = '#c7c7c7'
		cell_paths = self.ap.get_cells(variant=[variant_control])
		component_to_data = {component: [] for component in components}

		for i, sim_dir in enumerate(cell_paths):

			# Listeners used
			sim_out_dir = os.path.join(sim_dir, 'simOut')
			mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

			# Get data
			for component in components:
				mass = mass_reader.readColumn(component)
				component_to_data[component].append(mass[-1] / mass[[0]])

		x = -2
		median, yerr = get_bar_parameters(
			component_to_data['dnaMass'], per_generation=False)
		axes[0].bar(
			x, median, color=color_control, yerr=yerr, error_kw=error_kw)

		median, yerr = get_bar_parameters(
			component_to_data['rnaMass'], per_generation=False)
		axes[1].bar(
			x, median, color=color_control, yerr=yerr, error_kw=error_kw)

		median, yerr = get_bar_parameters(
			component_to_data['proteinMass'], per_generation=False)
		axes[2].bar(
			x, median, color=color_control, yerr=yerr, error_kw=error_kw)

		################################################################
		# Plot mass of experimental variant
		color_experiment = '#6b6ecf'
		color_highlight = '#b5cf6b' # #e7ba52
		seeds = self.ap.n_seed
		generations = self.ap.n_generation
		component_to_data = {component: np.nan * np.ones((seeds, generations))
			for component in components}

		for seed in range(seeds):
			for generation in range(generations):

				# Get cell
				cell_paths = self.ap.get_cells(
					variant=[variant_experiment],
					seed=[seed],
					generation=[generation],
					)

				if len(cell_paths) == 0:
					continue

				sim_dir = cell_paths[0]
				sim_out_dir = os.path.join(sim_dir, 'simOut')

				# Check for early termination
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

				time = main_reader.readColumn('time')
				mass = mass_reader.readColumn('dryMass')

				hit_3_hour_limit = time[-1] - time[0] < 3 * 60 * 60
				insufficient_mass_accumulation = mass[-1] / mass[0] < 1.1

				if hit_3_hour_limit and insufficient_mass_accumulation:
					print(f'skip: seed {seed}, generation {generation}')
					continue

				# Calculate mass fold change
				for component in components:
					mass = mass_reader.readColumn(component)

					if component == 'dnaMass' and mass[0] == 0:
						mass += 0.001
					component_to_data[component][seed, generation]\
						= mass[-1] / mass[0]

		x = np.arange(generations)
		median, yerr = get_bar_parameters(component_to_data['dnaMass'])
		color = [color_experiment for _ in range(generations)]
		color[np.where(median <= 1.1)[0][0]] = color_highlight
		axes[0].bar(x, median, color=color, yerr=yerr, error_kw=error_kw)

		median, yerr = get_bar_parameters(component_to_data['rnaMass'])
		color = [color_experiment for _ in range(generations)]
		color[np.where(median <= 1.1)[0][0]] = color_highlight
		axes[1].bar(x, median, color=color, yerr=yerr, error_kw=error_kw)

		median, yerr = get_bar_parameters(component_to_data['proteinMass'])
		color = [color_experiment for _ in range(generations)]
		color[np.where(median <= 1.1)[0][0]] = color_highlight
		axes[2].bar(x, median, color=color, yerr=yerr, error_kw=error_kw)

		################################################################
		# Format
		fontsize_heading = 11
		fontsize_label = 9
		fontsize_legend = 7

		# plt.suptitle('Accumulation of Mass of Cellular Components',
		# 	fontsize=fontsize_heading)
		axes[0].set_title('DNA', fontsize=fontsize_label)
		axes[1].set_title('RNA', fontsize=fontsize_label)
		axes[2].set_title('Protein', fontsize=fontsize_label)
		axes[2].set_xlabel('Generation', fontsize=fontsize_legend)

		# axes[1].set_ylabel('Mass Fold\nChange', fontsize=fontsize_legend)

		for ax in axes:
			# ax.set_ylabel('Mass Fold Change', fontsize=fontsize_label)
			ax.set_yticks([0, 1, 2])
			ax.set_xticks(np.arange(0, generations, 2))
			ax.set_xlim([-3, 20])
			ax.axvline(-1, color='k')
			ax.axhline(1, color='k', linestyle='--', linewidth=0.5)
			ax.tick_params(axis='both', which='major', labelsize=7)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		plt.tight_layout(h_pad=0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
