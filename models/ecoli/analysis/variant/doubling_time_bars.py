"""
Plots doubling times for variants in the tRNA synthetase kcat study.
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

		# Analyze variants 0 (control) and 6 (ArgRS)
		if not (0 in variants and 6 in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 6.')
			return

		variant_experiment = 6

		# Get doubling times
		data = {variant: [] for variant in variants}

		variant = 0
		cell_paths = self.ap.get_cells(variant=[variant])
		for sim_dir in cell_paths:
			main_reader = TableReader(os.path.join(sim_dir, 'simOut', 'Main'))
			time = main_reader.readColumn('time') / 60
			data[variant].append(time[-1] - time[0])

		variant = variant_experiment
		matrix = np.nan * np.ones((self.ap.n_seed, self.ap.n_generation))
		n_non_viable = 0
		for seed in range(self.ap.n_seed):
			for generation in range(self.ap.n_generation):
				sim_dirs = self.ap.get_cells(
					variant=[variant], seed=[seed], generation=[generation])

				if len(sim_dirs) > 0:
					sim_dir = sim_dirs[0]

					# Mass accumulation
					mass_reader = TableReader(
						os.path.join(sim_dir, 'simOut', 'Mass'))
					dry_mass = mass_reader.readColumn('dryMass')
					mass_fold_change = dry_mass[-1] / dry_mass[0]
					mass_doubled = mass_fold_change > 1.5

					# Simulated doubling time
					main_reader = TableReader(
						os.path.join(sim_dir, 'simOut', 'Main'))
					time = main_reader.readColumn('time') / 60
					tau = time[-1] - time[0]

					# Consider viable sims as having either: doubled
					# their cell mass, or hit the 180 min upper limit
					viable = mass_doubled or tau >= 179

					if viable:
						matrix[seed, generation] = tau
					else:
						n_non_viable += 1
		
		n_not_simulated = np.isnan(matrix).sum()
		n_total_sims = self.ap.n_seed * self.ap.n_generation
		print('Variant 6'
			f'\n{n_total_sims} total sims'
			f'\n\t{n_not_simulated} not simulated'
			f'\n\t\t{n_non_viable} non-viable sims'
			f'\n\t\t{n_not_simulated - n_non_viable} sims downstream'
			f'\n\t{n_total_sims - n_not_simulated} simulated'
			)

		data[variant] = matrix

		################################################################
		# Plot
		fig, axes = plt.subplots(2, 1, figsize=(3, 1.75),
			gridspec_kw={'height_ratios': [1, 2]}, sharex=True)
		ax_n, ax_bar = axes

		color_control = '#c7c7c7'
		color_experiment = '#6b6ecf'
		error_kw = {'lw': 0.5, 'capsize': 2, 'capthick': 0.5}

		# Control
		control_median = np.nanmedian(data[0])
		q1 = np.nanquantile(data[0], 0.25)
		q3 = np.nanquantile(data[0], 0.75)
		control_yerr = np.reshape([
			control_median - q1,
			q3 - control_median,
			], (2, 1))
		ax_bar.bar(-2, control_median, color=color_control, yerr=control_yerr, error_kw=error_kw)

		# Experiment
		n_sims = []
		for x in range(data[variant_experiment].shape[1]):
			data_col = data[variant_experiment][:, x]
			y = np.nanmedian(data_col)
			q1 = np.nanquantile(data_col, 0.25)
			q3 = np.nanquantile(data_col, 0.75)
			yerr = np.reshape([y - q1, q3 - y], (2, 1)) 
			ax_bar.bar(x, y, color=color_experiment, yerr=yerr, error_kw=error_kw)

			if x == 0:
				print(f'Gen 0 / Control: {y / control_median:.1f}')

			n_sims.append(np.logical_not(np.isnan(data_col)).sum())

		# Number of sims per bar
		ax_n.plot(
			np.arange(len(n_sims)),
			n_sims,
			color='k',
			marker='o',
			markeredgewidth=0,
			markersize=3,
			linewidth=0.5,
			)

		# Format
		# ax_n.set_title('Doubling Times', fontsize=7)
		ax_n.spines['bottom'].set_visible(False)
		# ax_n.set_ylabel('Number of Simulations', fontsize=7)
		ax_n.set_yticks([3, 10])
		ax_n.set_ylim([2, 11])
		ax_n.tick_params(axis='x', which='both', bottom=False)     

		ax_bar.set_xlabel('Generation', fontsize=7)
		# ax_bar.set_ylabel('Doubling Time (min)', fontsize=7)
		ax_bar.axvline(-1, color='k')
		ax_bar.set_xticks(np.arange(0, self.ap.n_generation, 2))
		ax_bar.set_yticks([0, 60, 120, 180])

		for ax in axes:
			ax.tick_params(axis='both', which='major', labelsize=7)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.set_xlim([-3, 20])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
