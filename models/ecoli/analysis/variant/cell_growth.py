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


GLC_ID = 'GLC[p]'
FLUX_UNITS = units.mmol / units.g / units.h
MASS_UNITS = units.fg
GROWTH_UNITS = MASS_UNITS / units.s


def plot_validation(mean, std, labels, val_rates, val_std, val_aa_ids, label, text):
	# Normalize simulation data by the control condition
	rate_mapping = {label: rate for label, rate in zip(labels, mean)}
	std_mapping = {label: std for label, std in zip(labels, std)}
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

	# Statistics
	r, p = pearsonr(val_rates, wcm_normalized_growth_rates)
	n = len(val_rates)

	plt.errorbar(val_rates, wcm_normalized_growth_rates,
		xerr=val_std, yerr=wcm_normalized_std, fmt='o', alpha=0.5,
		label=f'{label} r={r:.2f} (p={p:.2g}, n={n})')

	if text:
		for aa, x, y in zip(val_aa_ids, val_rates, wcm_normalized_growth_rates):
			plt.text(x, 0.01 + y, aa, ha='center', fontsize=6)

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

		aa_ids = sim_data.molecule_groups.amino_acids
		glc_mw = sim_data.getter.get_mass(GLC_ID)

		variant_lengths = []
		variant_counts = []
		variant_growth_rates = []
		variant_elong_rates = []
		variant_glc_yields = []
		labels = []
		for variant in variants:
			lengths = []
			growth_rates = []
			elong_rates = []
			glc_yields = []
			count = 0
			for sim_dir in ap.get_cells(variant=[variant]):
				glc_yield = 0
				try:
					sim_out_dir = os.path.join(sim_dir, "simOut")

					# Listeners used
					main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
					mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
					ribosome_reader = TableReader(os.path.join(sim_out_dir, 'RibosomeData'))
					fba_results = TableReader(os.path.join(sim_out_dir, 'FBAResults'))

					# Load data
					time = main_reader.readColumn('time')
					cycle_length = time[-1] - time[0]
					growth_rate = mass_reader.readColumn('instantaneous_growth_rate')[1:].mean()
					elong_rate = ribosome_reader.readColumn('effectiveElongationRate')[1:].mean()
					ex_molecules = fba_results.readAttribute('externalMoleculeIDs')
					if 'GLC[p]' in ex_molecules:
						glc_idx = ex_molecules.index(GLC_ID)
						glc_flux = FLUX_UNITS * -fba_results.readColumn('externalExchangeFluxes')[:, glc_idx]
						yields = 1 / glc_flux.asNumber()

						# The above is not a true yield but useful for comparison until uptake problems are fixed
						# The commented code below will be a better yield measure (g cell / g glc) in the future
						# growth = GROWTH_UNITS * mass_reader.readColumn('growth') / main_reader.readColumn('timeStepSec')
						# dry_mass = MASS_UNITS * mass_reader.readColumn('dryMass')
						# yields = units.strip_empty_units(growth / (glc_flux * glc_mw * dry_mass))

						glc_yield = yields[np.isfinite(yields)].mean()
				except Exception as e:
					print(f'Exception reading Main/time, Mass/instantaneous_growth_rate,'
						  f' RibosomeData/effectiveElongationRate, FBAResults/externalMoleculeIDs,'
						  f' or FBAResults/externalExchangeFluxes: {e!r}')
					cycle_length = 0
					growth_rate = 0
					elong_rate = 0

				# Filter out cell cycle lengths that are too short (likely failed)
				# TODO: better way to test for failure
				# TODO: also check for long cells that are going to fail
				if cycle_length / 60 < 15:
					lengths.append(np.inf)
				else:
					lengths.append(cycle_length / 60)
					count += 1
				growth_rates.append(growth_rate)
				elong_rates.append(elong_rate)
				glc_yields.append(glc_yield)

			variant_lengths.append(lengths)
			variant_counts.append(count)
			variant_growth_rates.append(growth_rates)
			variant_elong_rates.append(elong_rates)
			variant_glc_yields.append(glc_yields)
			labels.append(aa_ids[variant][:-3])

		all_lengths = np.vstack(variant_lengths)
		mean_lengths = np.array([np.mean(row[np.isfinite(row) & (row > 0)]) for row in all_lengths])
		all_growth_rates = np.vstack(variant_growth_rates) * 3600
		mean_growth_rates = all_growth_rates.mean(axis=1)
		std_growth_rates = all_growth_rates.std(axis=1)
		all_elong_rates = np.vstack(variant_elong_rates)
		mean_elong_rates = all_elong_rates.mean(axis=1)
		std_elong_rates = all_elong_rates.std(axis=1)
		all_glc_yields = np.vstack(variant_glc_yields)
		mean_glc_yields = all_glc_yields.mean(axis=1)
		std_glc_yields = all_glc_yields.std(axis=1)

		# Load validation growth rates
		all_aa_ids = {aa[:-3] for aa in aa_ids}
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

		# Create plots
		plt.figure(figsize=(16, 10))
		gs = gridspec.GridSpec(5, 3)

		## Bar plot of cell cycle lengths
		ax = plt.subplot(gs[0, 0])
		plt.bar(variants, mean_lengths)
		plt.ylabel('Average cell cycle length (min)', fontsize=8)
		remove_border(ax, bottom=True)

		## Bar plot of growth rates
		ax = plt.subplot(gs[1, 0])
		plt.bar(variants, mean_growth_rates, yerr=std_growth_rates)
		plt.ylabel('Average growth rate (1/hr)', fontsize=8)
		remove_border(ax, bottom=True)

		## Bar plot of elongation rates
		ax = plt.subplot(gs[2, 0])
		plt.bar(variants, mean_elong_rates, yerr=std_elong_rates)
		plt.ylabel('Average elongation rate (AA/s)', fontsize=8)
		remove_border(ax, bottom=True)

		## Bar plot of glucose yield
		ax = plt.subplot(gs[3, 0])
		plt.bar(variants, mean_glc_yields, yerr=std_glc_yields)
		plt.ylabel('Glc yield (1 / uptake)', fontsize=8)
		remove_border(ax, bottom=True)

		## Bar plot of valid simulations
		ax = plt.subplot(gs[4, 0])
		plt.bar(variants, variant_counts)
		plt.ylabel('Number of variants', fontsize=8)
		plt.xticks(variants, labels, rotation=45, fontsize=6, ha='right')
		remove_border(ax)

		## Validation comparison for each amino acid addition
		if metadata.get('variant', '') == 'add_one_aa':
			ax = plt.subplot(gs[:, 1:])

			# Plot datasets to compare against validation
			plot_validation(mean_growth_rates, std_growth_rates, labels,
				val_normalized_growth_rates, val_normalized_std, val_aa_ids,
				'Growth rate', text=True)
			plot_validation(mean_elong_rates, std_elong_rates, labels,
				val_normalized_growth_rates, val_normalized_std, val_aa_ids,
				'Elong rate', text=False)
			plot_validation(mean_glc_yields, std_glc_yields, labels,
				val_normalized_growth_rates, val_normalized_std, val_aa_ids,
				'Glc yield', text=False)

			# Plot formatting
			plt.legend(fontsize=8)
			remove_border(ax)
			ax.tick_params(axis='x', labelsize=6)
			plt.xlabel('Validation growth rate\n(Normalized to minimal media)')
			plt.ylabel('Simulation data\n(Normalized to minimal media)')

			# Plot y=x diagonal
			x_min, x_max = plt.xlim()
			y_min, y_max = plt.ylim()
			min_rate = min(x_min, y_min)
			max_rate = max(x_max, y_max)
			plt.plot([min_rate, max_rate], [min_rate, max_rate], '--k')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
