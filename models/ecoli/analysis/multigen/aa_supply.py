"""
Show factors that contribute to the supply rates for amino acids.  Traces for
'supply vs expected', 'normalized enzyme conc' and 'normalized AA conc' should
stay around 1 for normal cells.  Can be helpful in troubleshooting issues with
synthesis enzyme production or increases/drops in the amino acid concentration
and the effect that both of these can have on the supply compared to expected.
"""

import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
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

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		times = read_stacked_columns(cell_paths, 'Main', 'time') / 60
		time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec')
		dry_mass = read_stacked_columns(cell_paths, 'Mass', 'dryMass')
		counts_to_mol = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar')
		supply = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply')
		enzyme_counts = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes')
		aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_aa_conc').T
		supply_fraction = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_fraction').T

		# Calculate derived quantities
		normalized_supply = (supply / time_step / dry_mass / expected_supply).T
		enzyme_conc = (enzyme_counts * counts_to_mol).T
		enzyme_conc /= enzyme_conc[:, 1:2]
		aa_conc /= aa_conc[:, 1:2]

		# Plot data
		plt.figure(figsize=(16, 12))
		n_subplots = n_aas + 1  # 1 for legend
		rows = int(np.ceil(np.sqrt(n_subplots)))
		cols = int(np.ceil(n_subplots / rows))
		gs = gridspec.GridSpec(nrows=rows, ncols=cols)

		## Plot data for each amino acid
		for i, (supply, enzymes, aa_conc, fraction) in enumerate(zip(normalized_supply, enzyme_conc, aa_conc, supply_fraction)):
			row = i // cols
			col = i % cols
			ax = plt.subplot(gs[row, col])

			ax.plot(times, supply, label='Supply vs expected')
			ax.plot(times, enzymes, label='Normalized enzyme conc')
			ax.plot(times, aa_conc, label='Normalized AA conc')
			ax.plot(times, fraction, label='Fraction max supply')

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
