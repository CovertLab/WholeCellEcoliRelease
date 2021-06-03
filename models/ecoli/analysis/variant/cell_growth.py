"""
Compare cell cycle times and growth rates across variants.  Useful as validation
with the add_one_aa variant and can also be used to compare variants in
remove_one_aa variant.
"""

import pickle
import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


def remove_border(ax, bottom=False):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	if bottom:
		ax.spines['bottom'].set_visible(False)
		ax.set_xticks([])
	ax.tick_params(axis='y', labelsize=6)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		variant_lengths = []
		variant_counts = []
		variant_rates = []
		labels = []
		for variant in variants:
			lengths = []
			rates = []
			count = 0
			for sim_dir in ap.get_cells(variant=[variant]):
				try:
					sim_out_dir = os.path.join(sim_dir, "simOut")

					# Listeners used
					main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
					mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

					# Load data
					time = main_reader.readColumn('time')
					cycle_length = time[-1] - time[0]
					growth_rate = mass_reader.readColumn('instantaneous_growth_rate')[1:].mean()
				except Exception:
					cycle_length = 0
					growth_rate = 0

				# Filter out cell cycle lengths that are too short (likely failed)
				# TODO: better way to test for failure
				# TODO: also check for long cells that are going to fail
				if cycle_length / 60 < 15:
					lengths.append(np.inf)
				else:
					lengths.append(cycle_length / 60)
					count += 1
				rates.append(growth_rate)

			variant_lengths.append(lengths)
			variant_counts.append(count)
			variant_rates.append(rates)
			labels.append(sim_data.molecule_groups.amino_acids[variant][:-3])

		all_lengths = np.vstack(variant_lengths)
		mean_lengths = np.array([np.mean(row[np.isfinite(row) & (row > 0)]) for row in all_lengths])
		all_rates = np.vstack(variant_rates) * 3600
		mean_rates = all_rates.mean(axis=1)
		std_rates = all_rates.std(axis=1)

		# Load validation growth rates
		all_aa_ids = {aa[:-3] for aa in sim_data.molecule_groups.amino_acids}
		val_control = validation_data.amino_acid_growth_rates['minimal']['mean']
		val_aa_ids = []
		val_normalized_growth_rates = []
		val_normalized_std = []
		for media, rates in validation_data.amino_acid_growth_rates.items():
			aa_id = media.split('_')[-1]
			if aa_id in all_aa_ids:
				val_aa_ids.append(aa_id)
				val_normalized_growth_rates.append(units.strip_empty_units(rates['mean'] / val_control))
				val_normalized_std.append(units.strip_empty_units(rates['std'] / val_control))
		val_normalized_growth_rates = np.array(val_normalized_growth_rates)
		val_normalized_std = np.array(val_normalized_std)

		# Normalize simulation rates by the control condition
		rate_mapping = {label: rate for label, rate in zip(labels, mean_rates)}
		std_mapping = {label: std for label, std in zip(labels, std_rates)}
		control_label = 'L-SELENOCYSTEINE'  # control because SEL is already included for uptake in minimal media
		wcm_control = rate_mapping.get(control_label, 1)
		wcm_normalized_growth_rates = np.array([
			rate_mapping.get(aa, 0) / wcm_control
			for aa in val_aa_ids
			])
		wcm_normalized_std = np.array([
			std_mapping.get(aa, 0) / wcm_control
			for aa in val_aa_ids
			])

		# Create plots
		plt.figure(figsize=(16, 8))
		gs = gridspec.GridSpec(3, 2)

		## Bar plot of cell cycle lengths
		ax = plt.subplot(gs[0, 0])
		plt.bar(variants, mean_lengths)
		plt.ylabel('Average cell cycle length (min)', fontsize=8)
		remove_border(ax, bottom=True)

		## Bar plot of growth rates
		ax = plt.subplot(gs[1, 0])
		plt.bar(variants, mean_rates, yerr=std_rates)
		plt.ylabel('Average growth rate (1/hr)', fontsize=8)
		remove_border(ax, bottom=True)

		## Bar plot of valid simulations
		ax = plt.subplot(gs[2, 0])
		plt.bar(variants, variant_counts)
		plt.ylabel('Number of variants', fontsize=8)
		plt.xticks(variants, labels, rotation=45, fontsize=6, ha='right')
		remove_border(ax)

		## Validation comparison for each amino acid addition
		if metadata.get('variant', '') == 'add_one_aa':
			# Statistics
			r, p = pearsonr(val_normalized_growth_rates, wcm_normalized_growth_rates)
			n = len(val_normalized_growth_rates)

			ax = plt.subplot(gs[:, 1])
			min_rate = min(val_normalized_growth_rates.min(), wcm_normalized_growth_rates.min())
			max_rate = max(val_normalized_growth_rates.max(), wcm_normalized_growth_rates.max())

			plt.errorbar(val_normalized_growth_rates, wcm_normalized_growth_rates,
				xerr=val_normalized_std, yerr=wcm_normalized_std, fmt='o')
			plt.plot([min_rate, max_rate], [min_rate, max_rate], '--k')
			for aa, x, y in zip(val_aa_ids, val_normalized_growth_rates, wcm_normalized_growth_rates):
				plt.text(x, 0.01 + y, aa, ha='center', fontsize=6)

			remove_border(ax)
			ax.tick_params(axis='x', labelsize=6)
			plt.xlabel('Validation growth rate\n(Normalized to minimal media)')
			plt.ylabel('Simulation growth rate\n(Normalized to minimal media)')
			plt.title(f'Growth comparison to validation data\nr={r:.2f} (p={p:.2g}, n={n})', fontsize=8)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
