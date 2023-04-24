"""
Plots impact of HisRS kcat drop.
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
from wholecell.utils import units


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# Check sims
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

		# Initialize plot
		fig, axes = plt.subplots(4, 1, figsize=(4.5, 7))
		color = '#6b6ecf'
		color_highlight = '#b5cf6b'
		color_control = '#c7c7c7'
		fontsize_ticks = 7
		fontsize_labels = 9
		fontsize_heading = 11

		################################################################
		# Doubling time
		cell_paths = self.ap.get_cells(variant=[variant_experiment])
		doubling_times = []
		seeds = []
		for sim_dir in cell_paths:
			sim_out_dir = os.path.join(sim_dir, 'simOut')
			main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
			time = main_reader.readColumn('time') / 60
			doubling_times.append(time[-1] - time[0])
			seeds.append(int(sim_dir.split('/')[-3]))
		mean_doubling_time = np.mean(doubling_times)

		# Seed focus
		seed = seeds[np.argmax(doubling_times)]
		print(f'Seed: {seed}')

		cell_paths = self.ap.get_cells(
			variant=[variant_experiment], seed=[seed])
		doubling_times = []
		division_times = []
		half_gen_times = []

		for sim_dir in cell_paths:
			sim_out_dir = os.path.join(sim_dir, 'simOut')
			main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
			time = main_reader.readColumn('time') / 60
			tau = time[-1] - time[0]
			doubling_times.append(tau)
			division_times.append(time[-1])
			half_gen_times.append(time[0] + tau / 2)

		ax = axes[0]
		x = half_gen_times
		height = doubling_times - mean_doubling_time
		bars = ax.bar(x, height, bottom=mean_doubling_time, color=color,
			width=12)

		# Highlight
		sorted_indexes = np.argsort(height)[::-1]
		for i in range(4):
			bars[sorted_indexes[i]].set_color(color_highlight)

		ax.set_title('Doubling Time',#\nRelative to Population Average',
			fontsize=fontsize_heading)
		ax.set_xlabel('Generation', fontsize=fontsize_labels)
		ax.set_ylabel('Time (min)', fontsize=fontsize_labels)
		ax.set_xticks(x)
		ax.set_xticklabels(np.arange(len(x)))
		ax.axhline(mean_doubling_time, color=color, linestyle='--',
			linewidth=0.5)
		ax.set_ylim([45, 75])
		ax.set_yticks([45, round(mean_doubling_time), 75])

		################################################################
		def draw_cell_cycle_boundaries(ax, y):
			for division_time in division_times:
				ax.axvline(division_time, color='k', linestyle='--',
					linewidth=0.5)

			for i, x in enumerate(half_gen_times):
				ax.text(x, y, f'Gen {i}', fontsize=fontsize_ticks,
					va='bottom', ha='center')
			return
		################################################################
		# Ribosome elongation rate
		control_cell_paths = self.ap.get_cells(variant=[variant_control])
		control_ribosome_rate = read_stacked_columns(
			control_cell_paths, 'RibosomeData', 'effectiveElongationRate',
			remove_first=True).reshape(-1)
		mean_ribosome_rate = np.mean(control_ribosome_rate)
		print(f'Control ribosome elongation rate, mean:'
			f' {mean_ribosome_rate:.1f}')

		# Seed focus
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True) / 60
		ribosome_rate = read_stacked_columns(
			cell_paths, 'RibosomeData', 'effectiveElongationRate',
			remove_first=True).reshape(-1)
		deviation = 0.40
		mask = ribosome_rate < (mean_ribosome_rate * (1 - deviation))

		ax = axes[1]

		# # Control average
		# ax.axhline(mean_ribosome_rate, color=color_control)

		# Highlight
		ax.plot(time[mask], ribosome_rate[mask],
			color=color_highlight, linewidth=0, marker='o', markersize=4,
			label=f'{int(deviation * 100)}%+ deviation from'
			f' {mean_ribosome_rate:.1f}')

		# All
		ax.plot(time[1:], ribosome_rate[1:], color=color, linewidth=1)
		print('Ribosome elongation rate:')
		print(f'\tmin:\t{ribosome_rate[1:].min():.1f}')
		i = np.argmin(ribosome_rate[1:])
		print(f'\tat:\t{time[1:][i][0]:.1f} min')

		# Format
		ax.set_title('Ribosome Elongation', fontsize=fontsize_heading)
		ax.set_xlabel('Time (min)', fontsize=fontsize_labels)
		ax.set_ylabel('Rate (aa/s/ribosome)', fontsize=fontsize_labels)
		ax.legend(loc='best', fontsize=fontsize_ticks)
		ax.set_ylim([8.5, 19])
		ax.set_yticks([9, 18])
		draw_cell_cycle_boundaries(ax, 18)

		################################################################
		# Aminoacylation rate
		free_trna = 'hisR-tRNA[c]'
		charged_trna = 'charged-hisR-tRNA[c]'
		synthetase = 'HISS-CPLX[c]'

		rna_data = sim_data.process.transcription.rna_data
		trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		trna_index = trnas.index(free_trna)

		charging_events = read_stacked_columns(
			control_cell_paths, 'TrnaCharging', 'charging_events',
			remove_first=True)[:, trna_index]
		cell_volume = 1e-15 * units.L * read_stacked_columns(
			control_cell_paths, 'Mass', 'cellVolume', remove_first=True).reshape(-1)
		time_step = units.s * read_stacked_columns(
			control_cell_paths, 'Main', 'timeStepSec', remove_first=True).reshape(-1)
		v_charging = (1
			/ sim_data.constants.n_avogadro
			/ cell_volume
			/ time_step
			* charging_events
			).asNumber(units.umol / units.L / units.s)
		control_v_charging = np.mean(v_charging)

		charging_events = read_stacked_columns(
			cell_paths, 'TrnaCharging', 'charging_events', remove_first=True
			)[:, trna_index]
		cell_volume = 1e-15 * units.L * read_stacked_columns(
			cell_paths, 'Mass', 'cellVolume', remove_first=True).reshape(-1)
		time_step = units.s * read_stacked_columns(
			cell_paths, 'Main', 'timeStepSec', remove_first=True).reshape(-1)
		v_charging = (1
			/ sim_data.constants.n_avogadro
			/ cell_volume
			/ time_step
			* charging_events
			).asNumber(units.umol / units.L / units.s)

		ax = axes[2]

		# # Control average
		# ax.axhline(control_v_charging, color=color_control)

		# Highlight
		ax.plot(time[mask], v_charging[mask],
			color=color_highlight, linewidth=0, marker='o', markersize=4)

		# All
		ax.plot(time[1:], v_charging[1:], color=color, linewidth=1)
		print('Aminoacylation rate:')
		i = np.argmin(np.abs(time[1:] - 20))
		print(f'\tstart (first 20 min):\t{v_charging[1:][:i + 1].mean():.1f}')
		print(f'\tmin:\t{v_charging[1:].min():.1f}')
		i = np.argmin(v_charging[1:])
		print(f'\tat:\t{time[1:][i][0]:.1f} min')
		i = np.argmin(ribosome_rate[1:])
		print(f'\tanother low at:\t{time[1:][i][0]:.1f} min, {v_charging[1:][i]} uM/s')

		# Format
		ax.set_title('hisR tRNA Aminoacylation', fontsize=fontsize_heading)
		ax.set_xlabel('Time (min)', fontsize=fontsize_labels)
		ax.set_ylabel('Rate (uM/s)', fontsize=fontsize_labels)
		ax.set_ylim([4, 8.5])
		ax.set_yticks([4, 8])
		draw_cell_cycle_boundaries(ax, 8)

		################################################################
		# HisRS concentration
		(n_synthetases,) = read_stacked_bulk_molecules(
			cell_paths, ([synthetase]), remove_first=True)
		c_synthetases = (1
			/ sim_data.constants.n_avogadro
			/ cell_volume
			* n_synthetases
			).asNumber(units.umol / units.L)

		ax = axes[3]

		# Highlight
		ax.plot(time[mask], c_synthetases[mask],
			color=color_highlight, linewidth=0, marker='o', markersize=4)

		# All
		ax.plot(time[1:], c_synthetases[1:], color=color, linewidth=1)
		print('HisRS concentration, highlighted time:')
		print(f'\tmax:\t{c_synthetases[mask].max():.2f}')
		print(f'\tmin:\t{c_synthetases[mask].min():.2f}')

		full_range = c_synthetases.max() - c_synthetases.min()
		portion_of_range = (c_synthetases[mask].max()
			- c_synthetases[mask].min()) / full_range
		print(f'\trange:\tlowest {portion_of_range * 100:.1f} %')
		i = np.argmin(c_synthetases[1:])
		print(f'\tat:\t{time[1:][i][0]:.1f} min')

		# Format
		ax.set_title('HisRS Cellular Abundance', fontsize=fontsize_heading)
		ax.set_xlabel('Time (min)', fontsize=fontsize_labels)
		ax.set_ylabel('Concentration (uM)', fontsize=fontsize_labels)
		ax.set_ylim([0.25, 0.6])
		ax.set_yticks([0.25, 0.6])
		draw_cell_cycle_boundaries(ax, 0.57)

		################################################################
		for ax in axes:
			ax.margins(0.01)

		xmin = np.inf
		xmax = 0
		for ax in axes:
			_xmin, _xmax = ax.get_xlim()
			if _xmin < xmin:
				xmin = _xmin
			if _xmax > xmax:
				xmax = _xmax

		for ax in axes:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.set_xlim([xmin, xmax])
			ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

		plt.tight_layout(h_pad=1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
