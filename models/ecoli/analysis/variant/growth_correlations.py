"""
Show the correlation between cell properties related to growth and the growth
rate.  Calculates the moving average of the growth rate for multiple windows
to show the time scale that some of these properties might act on when affecting
growth.  Many properties are correlated with each other so high correlation
with growth does not imply a direct effect and higher or lower levels of growth
might also affect the cell properties.

TODO:
	show scatter or 2D density plot for any attributes vs growth
"""

import os
import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader


# +1 to make the window balanced on either side of the point of interest
MA_WINDOWS = np.array([10, 100, 1000]) + 1
MA_START = np.array([10, 50, 100, 500, 1000, 5000])
REPLACEMENT_NAMES = {'L-ASPARTATE': 'ASP', 'L-ALPHA-ALANINE': 'ALA'}
REMOVED_NAMES = {'L-SELENOCYSTEINE'}


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_bar(self, ax, cell_property, growth_rates, xlabel, tick_labels, growth_labels=None):
		n_subproperties = cell_property.shape[0]
		n_growth_rates = len(growth_rates)
		all_r = np.zeros((n_subproperties, n_growth_rates))
		for i, property_ in enumerate(cell_property):
			for j, growth in enumerate(growth_rates):
				filter_mask = np.isfinite(property_[:len(growth)]) & np.isfinite(growth)  # filter to prevent pearsonr ValueError
				if np.any(filter_mask):
					all_r[i, j] = stats.pearsonr(property_[:len(growth)][filter_mask], growth[filter_mask])[0]

		width = 0.8 / n_growth_rates
		offsets = np.arange(n_growth_rates) * width - 0.4 + width / 2
		x = np.arange(n_subproperties)
		for i, (r, offset) in enumerate(zip(all_r.T, offsets)):
			label = None if growth_labels is None else growth_labels[i]
			ax.bar(x + offset, r, width, label=label)

		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xticks(x)
		ax.set_xticklabels(tick_labels, fontsize=6, rotation=45, ha='right')
		ax.set_xlabel(xlabel, fontsize=8)

		return ax.get_legend_handles_labels()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = np.array([aa[:-3] for aa in sim_data.molecule_groups.amino_acids])
		for original, replacement in REPLACEMENT_NAMES.items():
			aa_ids[aa_ids == original] = replacement
		removed_mask = np.array([aa not in REMOVED_NAMES for aa in aa_ids], dtype=bool)

		variants = self.ap.get_variants()
		n_variants = len(variants)

		scaling = 3
		n_cols = 9
		plt.figure(figsize=(scaling*n_cols, scaling*n_variants))
		gs = gridspec.GridSpec(nrows=n_variants, ncols=n_cols)

		for row, variant in enumerate(variants):
			cell_paths = self.ap.get_cells(variant=[variant])

			# Load attributes
			unique_molecule_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
			unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
			rnap_idx = unique_molecule_ids.index('active_RNAP')
			ribosome_idx = unique_molecule_ids.index('active_ribosome')

			# Load data
			average_time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec',
				remove_first=True, ignore_exception=True).mean()
			sim_time = read_stacked_columns(cell_paths, 'Main', 'time',
				remove_first=True, ignore_exception=True).squeeze()
			counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True, ignore_exception=True).squeeze()
			growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
				remove_first=True, ignore_exception=True).squeeze() * 3600
			aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc',
				remove_first=True, ignore_exception=True).T
			uncharged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'uncharged_trna_conc',
				remove_first=True, ignore_exception=True).T
			charged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'charged_trna_conc',
				remove_first=True, ignore_exception=True).T
			synthetase_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'synthetase_conc',
				remove_first=True, ignore_exception=True).T
			aa_supply_enzymes_fwd = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_fwd',
				remove_first=True, ignore_exception=True).T * counts_to_molar
			aa_supply_enzymes_rev = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_rev',
				remove_first=True, ignore_exception=True).T * counts_to_molar
			ppgpp_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True, ignore_exception=True).squeeze()
			unique_mol_counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				remove_first=True, ignore_exception=True)

			# Derived values
			fraction_charged = charged_trna / (uncharged_trna + charged_trna)
			active_ribosome_conc = unique_mol_counts[:, ribosome_idx] * counts_to_molar
			active_rnap_conc = unique_mol_counts[:, rnap_idx] * counts_to_molar

			# Apply moving average to growth to calculate correlation between current value and future growth
			conv_arrays = []
			ma_labels = ['Current growth']
			for window in MA_WINDOWS:
				# Create spacing on bar plots
				conv_arrays.append(None)
				ma_labels.append('')

				for start in MA_START:
					if start + window < len(growth):
						ar = np.ones(start+window)
						ar[:start] = 0
						conv_arrays.append(ar)
						ma_labels.append(f'Start: +{start*average_time_step:.0f} s, window: {(window-1)*average_time_step:.0f} s')
			all_growth = [growth] + [
				np.convolve(growth, conv / conv.sum(), mode='valid') if conv is not None else np.array([])
				for conv in conv_arrays
				]

			# Create subplots for this variant
			## Plot combined data
			combined_data = np.vstack((ppgpp_conc, active_ribosome_conc, active_rnap_conc, sim_time))
			combined_labels = ['ppGpp', 'Ribosome', 'RNAP', 'Time']
			ax = plt.subplot(gs[row, 0])
			handles, labels = self.plot_bar(ax, combined_data, all_growth,
				'Other conc', combined_labels, growth_labels=ma_labels)
			ax.set_ylabel(f'Variant {variant}\nGrowth correlation', fontsize=8)

			## Show legend for each moving average window of growth rates
			ax = plt.subplot(gs[row, -1])
			ax.axis('off')
			ax.legend(handles, labels, loc='center', frameon=False, fontsize=5, ncol=2)

			## Plot cell properties that exist for all amino acids
			per_aa_plot_data = [
				(aa_conc, 'Amino acid conc'),
				(uncharged_trna, 'Uncharged tRNA conc'),
				(charged_trna, 'Charged tRNA conc'),
				(fraction_charged, 'Fraction charged'),
				(synthetase_conc, 'Synthetase conc'),
				(aa_supply_enzymes_fwd, 'Forward enzyme conc'),
				(aa_supply_enzymes_rev, 'Reverse enzyme conc'),
				]
			for col, (property_, label) in enumerate(per_aa_plot_data):
				self.plot_bar(plt.subplot(gs[row, col + 1]), property_[removed_mask, :], all_growth, label, aa_ids[removed_mask])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
