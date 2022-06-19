"""
Plot growth rate time series data for each lineage for all variants normalized
by the control variant.  This is expected to be run with add_one_aa or
remove_one_aa variants and compared with the scatter plot shown in cell_growth
variant analysis plot.
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


CONTROL_LABEL = 'L-SELENOCYSTEINE'  # control because SEL is already included for uptake in minimal media
MOVING_WINDOW = 301  # Needs to be an odd number for even contribution from both sides
HALF_WINDOW = MOVING_WINDOW // 2


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		variants = self.ap.get_variants()
		n_variants = len(variants)

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		aa_ids = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids]

		# Get average growth rate in the control simulation if it was run
		control_variant = aa_ids.index(CONTROL_LABEL)
		if control_variant in variants:
			cell_paths = self.ap.get_cells(variant=[control_variant])
			control_growth_rate = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True, ignore_exception=True).mean()
		else:
			control_growth_rate = 1

		# Calculate expected increase in normalized growth rate for each amino acid addition
		validation_control = validation_data.amino_acid_growth_rates['minimal']['mean']
		validation_references = {}
		for media, values in validation_data.amino_acid_growth_rates.items():
			aa_id = media.split('_')[-1]
			if aa_id == 'minimal':
				aa_id = CONTROL_LABEL
			validation_references[aa_id] = units.strip_empty_units(values['mean'] / validation_control)

		n_rows = int(np.ceil(np.sqrt(n_variants)))
		n_cols = int(np.ceil(n_variants / n_rows))

		plt.figure(figsize=(3*n_rows, 3*n_cols))
		gs = gridspec.GridSpec(n_rows, n_cols)
		for i, variant in enumerate(variants):
			aa_id = aa_ids[variant]
			row = i // n_rows
			col = i % n_rows
			ax = plt.subplot(gs[row, col])
			cmap = plt.get_cmap('tab10')

			# Plot moving average of growth rate and horizontal line of average for each initial seed
			for seed in range(self.ap.n_seed):
				color = cmap(seed)
				cell_paths = self.ap.get_cells(variant=[variant], seed=[seed])
				times = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True) / 3600
				growth_rates = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
					remove_first=True, ignore_exception=True).squeeze()
				padded_growth = np.hstack((
					np.full(HALF_WINDOW, growth_rates[:HALF_WINDOW].mean()),
					growth_rates,
					np.full(HALF_WINDOW, growth_rates[-HALF_WINDOW:].mean())
					))
				averaged_rates = np.convolve(padded_growth, np.ones(MOVING_WINDOW) / MOVING_WINDOW, mode='valid')
				ax.plot(times, averaged_rates / control_growth_rate, '-', color=color, alpha=0.4)
				ax.axhline(growth_rates.mean() / control_growth_rate, color=color, linestyle='--', linewidth=0.5, alpha=0.8)

			# Add reference lines
			ax.axhline(1, linestyle='--', color='k', linewidth=0.5, alpha=0.4)
			if control_variant in variants:
				if aa_id in validation_references:
					ax.axhline(validation_references[aa_id], linestyle='--', color='k', linewidth=0.5, alpha=0.8)
				ylabel = 'Normalized growth to control'
			else:
				ylabel = 'Growth rate'

			if col == 0:
				ax.set_ylabel(ylabel, fontsize=8)

			# Format axes for better display
			ax.tick_params(axis='both', labelsize=6)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			if row == n_rows - 1:
				ax.set_xlabel('Time (hr)', fontsize=8)
			else:
				ax.spines['bottom'].set_visible(False)
				ax.set_xticks([])
			ax.set_title('Control' if variant == control_variant else aa_id, fontsize=8)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
