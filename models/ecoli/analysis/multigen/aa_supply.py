"""
Show factors that contribute to the supply rates for amino acids.  Traces for
'supply vs expected', 'normalized enzyme conc' and 'normalized AA conc' should
stay around 1 for normal cells.  Can be helpful in troubleshooting issues with
synthesis enzyme production or increases/drops in the amino acid concentration
and the effect that both of these can have on the supply compared to expected.
For more detailed analysis of the enzyme expression, look at the
aa_supply_enzymes plot.
"""

import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Amino acid data
		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		# Expected supply data
		media = sim_data.conditions[sim_data.condition]['nutrients']
		expected_supply = (sim_data.translation_supply_rate[media] * sim_data.constants.n_avogadro).asNumber(1 / units.fg / units.s)

		cell_paths = self.ap.get_cells()

		# Load data
		times = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True) / 60
		time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec', remove_first=True)
		dry_mass = read_stacked_columns(cell_paths, 'Mass', 'dryMass', remove_first=True)
		counts_to_mol = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar', remove_first=True)
		supply = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply', remove_first=True)
		synthesis = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_synthesis', remove_first=True)
		imported = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_import', remove_first=True)
		aas_used = read_stacked_columns(cell_paths, 'GrowthLimits', 'aasUsed', remove_first=True)
		enzyme_counts_fwd = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_fwd', remove_first=True)
		enzyme_counts_rev = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_rev', remove_first=True)
		aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_aa_conc', remove_first=True).T
		supply_fraction_fwd = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_fraction_fwd', remove_first=True).T
		supply_fraction_rev = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_fraction_rev', remove_first=True).T

		# Calculate derived quantities
		normalized_supply = (supply / time_step / dry_mass / expected_supply).T
		normalized_synthesis = (synthesis / time_step / dry_mass / expected_supply).T
		normalized_imported = (imported / time_step / dry_mass / expected_supply).T
		normalized_use = (aas_used / time_step / dry_mass / expected_supply).T
		enzyme_conc_fwd = (enzyme_counts_fwd * counts_to_mol).T
		enzyme_conc_rev = (enzyme_counts_rev * counts_to_mol).T
		enzyme_conc_fwd /= enzyme_conc_fwd[:, 0:1]
		enzyme_conc_rev /= enzyme_conc_rev[:, 0:1]
		aa_conc /= aa_conc[:, 0:1]

		# Plot data
		plt.figure(figsize=(16, 12))
		n_subplots = n_aas + 1  # 1 for legend
		rows = int(np.ceil(np.sqrt(n_subplots)))
		cols = int(np.ceil(n_subplots / rows))
		gs = gridspec.GridSpec(nrows=rows, ncols=cols)

		## Plot data for each amino acid
		for i, (supply, enzymes_fwd, enzymes_rev, aa_conc, fraction_fwd, fraction_rev, synthesis, imported, use) in enumerate(zip(
				normalized_supply, enzyme_conc_fwd, enzyme_conc_rev, aa_conc, supply_fraction_fwd,
				supply_fraction_rev, normalized_synthesis, normalized_imported, normalized_use)):
			row = i // cols
			col = i % cols
			ax = plt.subplot(gs[row, col])

			ax.plot(times, supply, label='Supply vs expected', alpha=0.5)
			ax.plot(times, enzymes_fwd, label='Normalized enzyme conc (forward)', alpha=0.5)
			ax.plot(times, enzymes_rev, label='Normalized enzyme conc (reverse)', alpha=0.5)
			ax.plot(times, aa_conc, label='Normalized AA conc', alpha=0.5)
			ax.plot(times, fraction_fwd, label='Fraction max supply (forward)', alpha=0.5)
			ax.plot(times, fraction_rev, label='Fraction max supply (reverse)', alpha=0.5)
			ax.plot(times, synthesis, label='Synthesis vs expected', alpha=0.5)
			ax.plot(times, imported, label='Imported vs expected', alpha=0.5)
			ax.plot(times, use, label='Translation use vs expected', alpha=0.5)
			ax.axhline(0, linestyle='--', linewidth=0.5, color='k', alpha=0.3)
			ax.axhline(1, linestyle='--', linewidth=0.5, color='k', alpha=0.3)

			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.set_title(aa_ids[i], fontsize=8)

			if row == rows - 1:
				ax.set_xlabel('Time (min)')
			if col == 0:
				ax.set_ylabel('Normalized value')

		## Display legend
		handles, labels = ax.get_legend_handles_labels()
		ax = plt.subplot(gs[-1, -1], frameon=False)
		ax.axis('off')
		ax.legend(handles, labels, loc='center', frameon=False)

		## Save plot
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
