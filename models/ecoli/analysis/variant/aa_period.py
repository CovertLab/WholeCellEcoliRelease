"""
Analysis of amino acid concentration periods when allosteric inhibition is removed.
Useful with remove_aa_inhibition variant.

TODO:
 - plot all aa for each variant
 - highlight mutants in _lag and _normalized plots
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate, signal

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants import remove_aa_inhibition
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


def switch(data, input_drop=5, input_ma=60, deriv_ma=300, splrep=False):
	if input_drop > 1:
		n_points = len(data) // input_drop
		n_drop = len(data) % input_drop
		if n_drop == 0:
			index = slice(None)
		else:
			index = slice(-n_drop)
		data = data[index].reshape(n_points, -1).mean(1)

	if input_ma:
		input_ma = min(len(data) - 1, input_ma)
		data = np.convolve(data, np.ones(input_ma) / input_ma, 'valid')

	t = np.arange(len(data)) / 3600 * input_drop

	if splrep:
		spline = interpolate.splrep(t, data, k=5)
		deriv = interpolate.splev(t, spline, der=2)
	else:
		spline = interpolate.CubicSpline(t, data)
		deriv = spline.derivative(2)(t)

	if deriv_ma:
		deriv_ma = min(len(deriv) - 1, deriv_ma)
		deriv = np.convolve(deriv, np.ones(deriv_ma) / deriv_ma, 'same')

	sign = np.sign(deriv)

	return np.sum(sign[:-1] != sign[1:])

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		variants = self.ap.get_variants()
		n_variants = len(variants)

		baseline = None
		all_corr = {}
		all_periods = {}
		all_periods_deriv = {}
		all_periods_mean_adjusted = {}
		all_conc_mean = {}
		all_conc_std = {}
		autocorrelate = lambda x: signal.correlate(x, x, method='fft')
		for variant in variants:
			total_times = []
			var_corr = []
			var_conc = []
			periods = []
			periods_deriv = []
			for seed in self.ap.get_seeds(variant):
				cell_paths = self.ap.get_cells(variant=[variant], seed=[seed])

				if not np.all([self.ap.get_successful(cell) for cell in cell_paths]):
					continue

				# Load data
				sim_time = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True).squeeze()
				aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc', remove_first=True).T

				corr = np.apply_along_axis(autocorrelate, 1, aa_conc)
				var_corr.append(corr[:, corr.shape[1]//2:] / corr.mean(1).reshape(-1, 1))

				total_time = (sim_time[-1] - sim_time[0]) / 3600
				above = aa_conc > aa_conc.mean(1).reshape(-1, 1)
				period_switches = np.sum(above[:, :-1] != above[:, 1:], axis=1)
				period = total_time / (period_switches / 2)
				period[period == np.inf] = np.nan
				total_times.append(total_time)
				var_conc.append(aa_conc)
				periods.append(period)

				deriv_switches = np.apply_along_axis(switch, 1, aa_conc)
				period_deriv = total_time / (deriv_switches / 2)
				period_deriv[period_deriv == np.inf] = np.nan
				periods_deriv.append(period_deriv)

			all_corr[variant] = var_corr

			stacked_conc = np.hstack(var_conc)
			all_conc_mean[variant] = stacked_conc.mean(1)
			all_conc_std[variant] = stacked_conc.std(1)
			all_periods[variant] = np.nanmean(np.vstack(periods), 0)
			all_periods_deriv[variant] = np.nanmean(np.vstack(periods_deriv), 0)

			mean_periods = []
			for aa_conc, total_time in zip(var_conc, total_times):
				above = aa_conc > all_conc_mean[variant].reshape(-1, 1)
				mean_switches = np.sum(above[:, :-1] != above[:, 1:], axis=1)
				period = total_time / (mean_switches / 2)
				period[period == np.inf] = np.nan
				mean_periods.append(period)
			all_periods_mean_adjusted[variant] = np.nanmean(np.vstack(mean_periods), 0)

		# Save average data for comparison across runs
		with open(f'{os.path.join(plotOutDir, plotOutFileName)}.tsv', 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			headers = ['Variant']
			for aa in aa_ids:
				headers += [f'{aa} mean', f'{aa} std', f'{aa} CV', f'{aa} period']
			writer.writerow(headers)

			for variant in variants:
				cols = [variant]
				for mean, std, period in zip(all_conc_mean[variant], all_conc_std[variant], all_periods[variant]):
					cols += [mean, std, std / mean, period]

				writer.writerow(cols)

		control_idx = 0
		all_conc_cv = {
			variant: all_conc_std[variant] / all_conc_mean[variant]
			for variant in variants
			}
		width = 0.4

		# Plot for each amino acid for all variants
		def plot_bar(data, ylabel, ylim=None, file_label=''):
			plt.figure(figsize=(1.8, 1.5))

			controls = [
				data.get(control_idx, np.zeros(n_aas))[aa_ids.index(remove_aa_inhibition.get_aa_and_ki_factor(variant)[0])]
				for variant in variants
				if variant != control_idx
				]
			mutants = [
				data[variant][aa_ids.index(remove_aa_inhibition.get_aa_and_ki_factor(variant)[0])]
				for variant in variants
				if variant != control_idx
				]
			xlabels = [
				remove_aa_inhibition.get_aa_and_ki_factor(variant)[0][:-3]
				for variant in variants
				if variant != control_idx
				]
			x = range(len(mutants))
			# TODO: plot value for all mutants not just specific to aa
			plt.bar(x, controls, width=-width, align='edge')
			plt.bar(x, mutants, width=width, align='edge')

			plt.xticks(x, xlabels, fontsize=8, rotation=45)
			plt.ylabel(ylabel, fontsize=8)
			if ylim:
				plt.ylim(ylim)

			plt.gca().tick_params(labelsize=8)

			self.remove_border()
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + file_label, metadata)
			plt.close('all')

		plot_bar(all_periods, 'Period')
		plot_bar(all_conc_cv, 'Coefficient of variation', file_label='_cv')
		plot_bar(all_periods_deriv, 'Period', ylim=(0, 1), file_label='_deriv')  # TODO: adjust/remove ylim if needed, hardcoded for paper
		plot_bar(all_periods_mean_adjusted, 'Period', file_label='_mean')

		mean_corr = {}
		for variant, var_corr in all_corr.items():
			max_n = max([corr.shape[1] for corr in var_corr])
			n = np.zeros(max_n)
			total_corr = np.zeros((n_aas, max_n))
			for corr in var_corr:
				n_timepoints = corr.shape[1]
				n[:n_timepoints] += 1
				total_corr[:, :n_timepoints] += corr
			mean_corr[variant] = total_corr / n

		control_corr = mean_corr.get(0)

		def plot(label='', normalized=False):
			_, axes = plt.subplots(nrows=n_variants, ncols=n_aas, figsize=(30, 2*n_variants))
			if n_variants == 1:
				axes = axes.reshape(1, -1)

			for row, all_corr, in enumerate(mean_corr.values()):
				for col, corr in enumerate(all_corr):
					ax = axes[row, col]

					if normalized and control_corr is not None:
						control_pad = 0.95  # approaches 0 at the longest lag which can cause ratio to explode
						n_overlap = int(min(len(corr), control_corr.shape[1]*control_pad))
						ax.plot(np.arange(n_overlap) / 3600, corr[:n_overlap] / control_corr[col, :n_overlap])
						ax.axhline(1, color='k', alpha=0.5, linestyle='--', linewidth=1)
					else:
						ax.plot(np.arange(len(corr)) / 3600, corr)
						if control_corr is not None:
							ax.plot(np.arange(control_corr.shape[1]) / 3600, control_corr[col, :])

					ax.tick_params(labelsize=6)

					if row == 0:
						ax.set_title(aa_ids[col], fontsize=6)

					if col == 0:
						ax.set_ylabel(f'Variant {variants[row]}\nautocorrelation', fontsize=6)

					if row == n_variants - 1:
						ax.set_xlabel('Lag (hr)', fontsize=6)

			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + label, metadata)
			plt.close('all')

		plot(label='_lag')
		plot(label='_normalized', normalized=True)


if __name__ == "__main__":
	Plot().cli()
