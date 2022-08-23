"""
Compare the growth rates in variants with varying synthesis parameters for each
amino acid to see how sensitive growth rate is to a range of parameters.
Useful with the aa_synthesis_sensitivity variant.
"""

import pickle

from matplotlib import colors, gridspec, pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants import aa_synthesis_sensitivity
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


# These are set in the variant and will need to be updated if there are changes to the media
GLT_INDEX = 5
CONTROL_INDEX = 19
PARAM_CONTROL_LABEL = 'L-SELENOCYSTEINE aa_kcats_fwd'  # This parameter should not affect growth rate


def calculate_sensitivity(data, variant, factors, attr, default=None):
	params = []
	values = []
	slopes = []
	if variant in data:
		for param, factor_results in data[variant].items():
			attrs = []
			for factor in factors:
				if factor == 1 and default is not None:
					value = default
				elif factor not in factor_results or not np.isfinite(value := factor_results[factor][attr]):
					break
				attrs.append(value)
			else:
				params.append(param)
				values.append(np.array(attrs))
				slopes.append(stats.linregress(np.log10(factors), attrs).slope)

	return np.array(params), np.array(values), np.array(slopes)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		variants = self.ap.get_variants()

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = sim_data.molecule_groups.amino_acids
		n_params = len(aa_synthesis_sensitivity.PARAMETERS)

		# Load simulation growth rates
		data = {}
		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		for variant in variants:
			media_index = aa_synthesis_sensitivity.get_media_index(variant, sim_data)
			aa_index = aa_synthesis_sensitivity.get_aa_index(variant, sim_data)
			param, factor = aa_synthesis_sensitivity.get_adjustment(variant)
			aa_adjusted = aa_ids[aa_index][:-3]
			param_label = f'{aa_adjusted} {param}'

			# Load data
			cells = self.ap.get_cells(variant=[variant])
			time_step = read_stacked_columns(cells, 'Main', 'timeStepSec',
				remove_first=True, ignore_exception=True).squeeze()
			growth_rate = read_stacked_columns(cells, 'Mass', 'instantaneous_growth_rate',
				remove_first=True, ignore_exception=True).mean() * 3600
			elong_rate = read_stacked_columns(cells, 'RibosomeData', 'effectiveElongationRate',
				remove_first=True, ignore_exception=True).mean()
			protein_growth = (read_stacked_columns(cells, 'Mass', 'proteinMass',
				fun=growth_function, ignore_exception=True).squeeze() / time_step).mean() * 3600
			rna_growth = (read_stacked_columns(cells, 'Mass', 'rnaMass',
				fun=growth_function, ignore_exception=True).squeeze() / time_step).mean() * 3600

			variant_data = {
				'Growth rate': growth_rate,
				'Elongation rate': elong_rate,
				'Protein growth rate': protein_growth,
				'RNA growth rate': rna_growth,
				}
			media_data = data.get(media_index, {})
			param_data = media_data.get(param_label, {})
			param_data[factor] = variant_data
			media_data[param_label] = param_data
			data[media_index] = media_data

		# Plot data
		plt.figure(figsize=(32, 32))
		gs = gridspec.GridSpec(3, 2)

		def subplot(label, data_index, subplot_index, ylabel=None, normalized=False):
			if data_index not in data:
				return

			plt.subplot(gs[subplot_index])
			selected_data = data[data_index]
			param_labels = list(selected_data.keys())
			label_indices = {param: i for i, param in enumerate(param_labels)}

			x = []
			y = []
			scale = []
			for param_label, factors in selected_data.items():
				for factor, values in factors.items():
					x.append(label_indices[param_label])
					denom = data.get(CONTROL_INDEX, {}).get(param_label, {}).get(factor, {}).get(label, 1) if normalized else 1
					y.append(values[label] / denom)
					scale.append(factor)
			scale = np.array(scale)
			edges = ['silver' if factor == 0 else 'black' for factor in scale]
			scale[scale == 0] = 1
			scatter = plt.scatter(x, y, c=scale, edgecolors=edges, cmap='RdBu', alpha=0.8, norm=colors.LogNorm())
			plt.colorbar(scatter)

			# Plot formatting
			ylabel = f'\n{ylabel}' if ylabel else ''
			normalized_text = '\n(Normalized to minimal media)' if normalized else ''
			self.remove_border()
			plt.xticks(range(len(param_labels)), param_labels, rotation=45, fontsize=6, ha='right')
			plt.yticks(fontsize=6)
			plt.ylabel(f'{label}{ylabel}{normalized_text}', fontsize=6)

			# Divide parameters by amino acid
			for vertical in range(n_params, np.max(x), n_params):
				plt.axvline(vertical - 0.5, color='gray', linestyle='--', linewidth=1, alpha=0.5)

		subplot('Growth rate', GLT_INDEX, (0, 0), normalized=True)
		plt.title('Effect of scaling synthesis parameters by a factor\n'
			'Color represents parameter scale factor\n'
			'Light gray circles have a value of 0 for the parameter', fontsize=6)
		subplot('Growth rate', GLT_INDEX, (1, 0), ylabel='Glt added')
		subplot('Growth rate', CONTROL_INDEX, (2, 0), ylabel='Minimal glc')
		subplot('Elongation rate', GLT_INDEX, (0, 1), normalized=True)
		subplot('Elongation rate', GLT_INDEX, (1, 1), ylabel='Glt added')
		subplot('Elongation rate', CONTROL_INDEX, (2, 1), ylabel='Minimal glc')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		if CONTROL_INDEX in data:
			# Identify groups of factors to split into subplots
			nonzero_factors = [f for f in aa_synthesis_sensitivity.FACTORS if f != 0]
			nonzero_factors_with_control = sorted(nonzero_factors + [1])
			increase_factors = [f for f in nonzero_factors_with_control if f >= 1]
			decrease_factors = [f for f in nonzero_factors_with_control if f <= 1]

			# Calculate the control growth rate
			params, all_rates, _ = calculate_sensitivity(data, CONTROL_INDEX, nonzero_factors, 'Growth rate')
			control_growth = all_rates[params == PARAM_CONTROL_LABEL, :]
			control_growth_rate = control_growth.mean()
			if not np.allclose(control_growth, control_growth_rate):
				raise ValueError('Control parameter results in variable growth rates.'
					' Run sims with no modified parameter or consider the mean.')

			# Calculate changes from baseline and max growth rates for each parameter
			highest_param_rates = all_rates[:, -1]
			lowest_param_rates = all_rates[:, 0]
			diff = highest_param_rates - lowest_param_rates
			highest_param_change = highest_param_rates - control_growth_rate
			lowest_param_change = lowest_param_rates - control_growth_rate
			max_rates = np.max(all_rates, 1)

			def slopes_plot(variant, factors, attr, axes, n_labeled=5, stds=1., control=None):
				params, all_rates, slopes = calculate_sensitivity(data, variant, factors, attr, default=control)

				slope_sort_idx = np.argsort(slopes)
				mean = slopes.mean()
				std = slopes.std()
				upper_limit = max(mean + std * stds, slopes[slope_sort_idx[-n_labeled-1]])
				lower_limit = min(mean - std * stds, slopes[slope_sort_idx[n_labeled]])

				bar_ax, trace_ax = axes

				# TODO: label params on x
				bar_ax.bar(range(len(slopes)), slopes[slope_sort_idx])
				bar_ax.set_ylabel(f'Slope of {attr} vs change in param', fontsize=8)
				bar_ax.tick_params(labelsize=8)
				self.remove_border(bar_ax)

				for param, rates, slope in zip(params, all_rates, slopes):
					if slope > upper_limit or slope < lower_limit:
						style = dict(alpha=0.6, linewidth=1, markersize=2, label=param)
					else:
						style = dict(color='k', alpha=0.3, linewidth=0.5, markersize=1)
					trace_ax.plot(np.log10(factors), rates, 'o-', **style)
				trace_ax.set_xlabel('log10 change in param', fontsize=8)
				trace_ax.set_ylabel(f'{attr}', fontsize=8)
				trace_ax.legend(fontsize=6, frameon=False)
				trace_ax.tick_params(labelsize=8)
				self.remove_border(trace_ax)

			_, axes = plt.subplots(2, 6, figsize=(30, 10))

			# Slopes of growth vs change in parameter
			# TODO: control growth rates for RNA and protein instead of overall growth rate
			slopes_plot(CONTROL_INDEX, nonzero_factors_with_control, 'Growth rate', axes[:, 0], control=control_growth_rate)
			slopes_plot(CONTROL_INDEX, nonzero_factors_with_control, 'Protein growth rate', axes[:, 1], control=control_growth_rate)
			slopes_plot(CONTROL_INDEX, nonzero_factors_with_control, 'RNA growth rate', axes[:, 2], control=control_growth_rate)
			slopes_plot(CONTROL_INDEX, increase_factors, 'Growth rate', axes[:, 3], control=control_growth_rate)
			slopes_plot(CONTROL_INDEX, decrease_factors, 'Growth rate', axes[:, 4], control=control_growth_rate)

			# Greatest changes from baseline in positive and negative directions
			ax = axes[0, 5]
			sort_idx = np.argsort(diff)
			x = np.arange(len(diff))
			ax.bar(x, highest_param_change[sort_idx])
			ax.bar(x, lowest_param_change[sort_idx])
			ax.set_ylabel('Change in growth from baseline (1/hr)', fontsize=8)
			ax.tick_params(labelsize=8)
			self.remove_border(ax)

			# Highest growth rates possible per param change
			ax = axes[1, 5]
			sort_idx = np.argsort(max_rates)
			x = np.arange(len(max_rates))
			ax.bar(x, max_rates[sort_idx])
			ax.set_ylabel('Highest growth rate per parameter (1/hr)', fontsize=8)
			ax.tick_params(labelsize=8)
			self.remove_border(ax)

			plt.tight_layout()
			exportFigure(plt, plotOutDir, f'{plotOutFileName}_aggregated', metadata)
			plt.close('all')
		else:
			print('Need to run additional variants for the aggregated plot.')


if __name__ == "__main__":
	Plot().cli()
