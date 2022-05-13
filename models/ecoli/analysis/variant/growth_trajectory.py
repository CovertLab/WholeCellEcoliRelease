"""
Show the trajectory of growth rates and cell composition measurements.  Useful
for variants with a shift.

TODO:
	- shares a lot of code with the multigen plot with the same name (make into a function?)
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, read_stacked_bulk_molecules
from wholecell.utils import units


def mean_std(data):
	data = np.hstack(data)
	return np.mean(data), np.std(data)

def plot(ax, x, y, sim_time=None, timeline=None, ma_time=None, xlabel=None, ylabel=None,
		label=None, background=False, markersize=3):
	# Markers for any timeline shifts
	if sim_time is not None and timeline is not None:
		for t, media in timeline:
			times_before_shift = (sim_time < t)[:len(x)]
			times_before_shift[0] = True
			marker_x = x[times_before_shift][-1:]
			marker_y = y[times_before_shift][-1:]
			ax.plot(marker_x, marker_y, 'rx')
			ax.text(marker_x, marker_y, media, fontsize=6)

	# Settings for different plots and start/end markers for main plot
	if background:
		kwargs = {'alpha': 0.2, 'linewidth': 0.6, 'color': 'black'}
	else:
		kwargs = {'alpha': 0.7, 'label': label}
		ax.plot(x[0], y[0], 'og', markersize=markersize)
		ax.plot(x[-1], y[-1], 'or', markersize=markersize)

	# Plot trace and time markers
	trace, = ax.plot(x, y, **kwargs)
	if ma_time is not None:
		time_hr = np.floor(ma_time / 3600)
		hour_markers = np.where(np.diff(time_hr))[0] + 1
		ax.plot(x[hour_markers], y[hour_markers], 'o', alpha=kwargs['alpha'],
			markeredgewidth=0, markersize=markersize*2, color=trace.get_color())

	# Format axes
	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel(ylabel, fontsize=8)
	ax.tick_params(labelsize=6)

def set_lim(ax, xmin=0, xmax=0.6, ymin=0, ymax=2):
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		variants = self.ap.get_variants()

		# Create plot
		_, main_axes = plt.subplots(3, 2, figsize=(8, 12))
		main_fig = plt.gcf().number
		_, trimmed_axes = plt.subplots(3, 2, figsize=(8, 12))
		trimmed_fig = plt.gcf().number


		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		average_data = {}
		for variant in variants:
			with open(self.ap.get_variant_kb(variant), 'rb') as f:
				sim_data = pickle.load(f)
			if sim_data.external_state.current_timeline_id:
				timeline = sim_data.external_state.saved_timelines[sim_data.external_state.current_timeline_id]
			else:
				timeline = []

			aa_ids = sim_data.molecule_groups.amino_acids
			aa_mws = sim_data.getter.get_masses(aa_ids).asNumber(units.fg / units.count)

			all_mass_means = []
			all_growth_means = []
			all_ratio_means = []
			all_ratio_ma = []
			all_rna_to_aa_ma = []
			all_growth_ma = []
			all_protein_growth_ma = []
			all_rna_growth_ma = []
			all_small_mol_growth_ma = []
			all_times = []
			all_times_ma = []
			all_growth = []
			all_ratio = []
			data_to_plot = False
			for seed in self.ap.get_seeds(variant):
				cell_paths = self.ap.get_cells(variant=[variant], seed=[seed])
				if len(cell_paths) == 0 or not self.ap.get_successful(cell_paths[-1]):
					continue

				# Load data
				sim_time = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True).squeeze()
				time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec', remove_first=True).squeeze() / 3600  # hr
				growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True).squeeze() * 3600  # 1/hr
				protein = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
				rna = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
				protein_growth = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', fun=growth_function).squeeze() / time_step
				rna_growth = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', fun=growth_function).squeeze() / time_step
				small_mol_growth = read_stacked_columns(cell_paths, 'Mass', 'smallMoleculeMass', fun=growth_function).squeeze() / time_step
				growth_means = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True, fun=np.mean).reshape(-1) * 3600
				protein_means = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True, fun=np.mean).reshape(-1)
				rna_means = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True, fun=np.mean).reshape(-1)
				mass_means = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True, fun=np.mean).reshape(-1)
				aas,  = read_stacked_bulk_molecules(cell_paths, (aa_ids,), remove_first=True)

				if len(np.unique(time_step)) > 1:
					raise ValueError('Check plot implementation to handle variable time step across sims.')

				moving_window = min(201, len(growth))
				convolution_array = np.ones(moving_window) / moving_window

				# Process data
				aa_mass = aas @ aa_mws
				ratio = rna / protein
				rna_to_aa = rna / (protein + aa_mass)
				ratio_means = rna_means / protein_means
				growth_ma = np.convolve(growth, convolution_array, mode='valid')
				ratio_ma = np.convolve(ratio, convolution_array, mode='valid')
				rna_to_aa_ma = np.convolve(rna_to_aa, convolution_array, mode='valid')
				protein_growth_ma = np.convolve(protein_growth, convolution_array, mode='valid')
				rna_growth_ma = np.convolve(rna_growth, convolution_array, mode='valid')
				small_mol_growth_ma = np.convolve(small_mol_growth, convolution_array, mode='valid')
				time_ma = np.convolve(sim_time, convolution_array, mode='valid')

				plot(main_axes[0, 0], mass_means, growth_means, background=True)
				plot(main_axes[1, 0], ratio_means, growth_means, background=True)
				plot(main_axes[2, 0], ratio_ma, growth_ma, background=True)
				plot(main_axes[0, 1], ratio_ma, protein_growth_ma, background=True)
				plot(main_axes[1, 1], ratio_ma, rna_growth_ma, background=True)
				plot(main_axes[2, 1], ratio_ma, small_mol_growth_ma, background=True)

				all_mass_means.append(mass_means)
				all_growth_means.append(growth_means)
				all_ratio_means.append(ratio_means)
				all_ratio_ma.append(ratio_ma)
				all_rna_to_aa_ma.append(rna_to_aa_ma)
				all_growth_ma.append(growth_ma)
				all_protein_growth_ma.append(protein_growth_ma)
				all_rna_growth_ma.append(rna_growth_ma)
				all_small_mol_growth_ma.append(small_mol_growth_ma)
				all_times.append(sim_time[moving_window-1:])  # offset to get first contribution of time point in ma for timeline shift
				all_times_ma.append(time_ma)
				all_growth.append(growth)
				all_ratio.append(ratio)
				data_to_plot = True

			if not data_to_plot:
				continue

			min_length = min([len(data) for data in all_growth_means])
			stacked_mass_means = np.vstack([data[:min_length] for data in all_mass_means]).mean(0)
			stacked_growth_means = np.vstack([data[:min_length] for data in all_growth_means]).mean(0)
			stacked_ratio_means = np.vstack([data[:min_length] for data in all_ratio_means]).mean(0)

			min_length_ma = min([len(data) for data in all_growth_ma])
			stacked_ratio_ma = np.vstack([data[:min_length_ma] for data in all_ratio_ma]).mean(0)
			stacked_growth_ma = np.vstack([data[:min_length_ma] for data in all_growth_ma]).mean(0)
			stacked_protein_growth_ma = np.vstack([data[:min_length_ma] for data in all_protein_growth_ma]).mean(0)
			stacked_rna_growth_ma = np.vstack([data[:min_length_ma] for data in all_rna_growth_ma]).mean(0)
			stacked_small_mol_growth_ma = np.vstack([data[:min_length_ma] for data in all_small_mol_growth_ma]).mean(0)
			stacked_times = np.vstack([data[:min_length_ma] for data in all_times]).mean(0)
			stacked_times_ma = np.vstack([data[:min_length_ma] for data in all_times_ma]).mean(0)

			for axes, tl in [(main_axes, timeline), (trimmed_axes, None)]:
				plot(axes[0, 0], stacked_mass_means, stacked_growth_means,
					xlabel='Average cell cycle mass (fg)', ylabel='Average cell cycle growth rate (1/hr)', label=variant)
				plot(axes[1, 0], stacked_ratio_means, stacked_growth_means,
					xlabel='Average cell cycle RNA/protein', ylabel='Average cell cycle growth rate (1/hr)', label=variant)
				plot(axes[2, 0], stacked_ratio_ma, stacked_growth_ma,
					ma_time=stacked_times_ma, sim_time=stacked_times, timeline=tl,
					xlabel='RNA/protein', ylabel='Growth rate (1/hr)', label=variant)
				plot(axes[0, 1], stacked_ratio_ma, stacked_protein_growth_ma,
					ma_time=stacked_times_ma, sim_time=stacked_times, timeline=tl,
					xlabel='RNA/protein', ylabel='Protein growth rate (1/hr)', label=variant)
				plot(axes[1, 1], stacked_ratio_ma, stacked_rna_growth_ma,
					ma_time=stacked_times_ma, sim_time=stacked_times, timeline=tl,
					xlabel='RNA/protein', ylabel='RNA growth rate (1/hr)', label=variant)
				plot(axes[2, 1], stacked_ratio_ma, stacked_small_mol_growth_ma,
					ma_time=stacked_times_ma, sim_time=stacked_times, timeline=tl,
					xlabel='RNA/protein', ylabel='Small molecule growth rate (1/hr)', label=variant)

			# Save average/std for output to a tsv
			average_data[variant] = {
				'Growth': mean_std(all_growth_ma),
				'R/P ratio': mean_std(all_ratio_ma),
				'R/(P+A) ratio': mean_std(all_rna_to_aa_ma),
				}

		for axes in [main_axes, trimmed_axes]:
			for ax in axes.reshape(-1):
				self.remove_border(ax)
				if len(variants) > 1:
					ax.legend(fontsize=6)

		plt.figure(main_fig)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.figure(trimmed_fig)
		set_lim(trimmed_axes[0, 0], xmin=0, xmax=4000)
		set_lim(trimmed_axes[1, 0])
		set_lim(trimmed_axes[2, 0], ymin=-0.5, ymax=3.5)
		set_lim(trimmed_axes[0, 1])
		set_lim(trimmed_axes[1, 1])
		set_lim(trimmed_axes[2, 1], ymin=-1, ymax=4.5)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_trimmed', metadata)

		plt.close('all')

		# Save average data for comparison across runs
		with open(f'{os.path.join(plotOutDir, plotOutFileName)}.tsv', 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			headers = ['Variant']
			for col in next(iter(average_data.values())):
				headers += [f'{col} mean', f'{col} std']
			writer.writerow(headers)

			for variant, data in average_data.items():
				cols = [variant]
				for col in data.values():
					cols += col

				writer.writerow(cols)


if __name__ == "__main__":
	Plot().cli()
