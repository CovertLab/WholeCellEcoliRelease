"""
Show the concentration of amino acid synthesis enzymes compared to the expected
values from the parca for basal and with amino acid conditions. Can be helpful
in troubleshooting issues with synthesis enzyme production and could be paired
with the aa_supply analysis plot.
"""

import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import CONC_UNITS
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


PLOT_UNITS = units.umol / units.L


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Amino acid data
		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		# Expected expression from parameter calculations
		metabolism = sim_data.process.metabolism
		enzyme_to_amino_acid = metabolism.enzyme_to_amino_acid
		with_aa_reference = metabolism.aa_supply_enzyme_conc_with_aa.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid
		basal_reference = metabolism.aa_supply_enzyme_conc_basal.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		times = read_stacked_columns(cell_paths, 'Main', 'time') / 60
		counts_to_mol = (CONC_UNITS * read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar')).asNumber(PLOT_UNITS)
		enzyme_counts = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes')

		# Calculate derived quantities
		start = times[0, 0]
		end = times[-1, 0]
		enzyme_conc = (enzyme_counts * counts_to_mol).T

		# Plot data
		plt.figure(figsize=(16, 12))
		n_subplots = n_aas + 1  # 1 for legend
		rows = int(np.ceil(np.sqrt(n_subplots)))
		cols = int(np.ceil(n_subplots / rows))
		gs = gridspec.GridSpec(nrows=rows, ncols=cols)

		## Plot data for each amino acid
		for i, (with_aa, basal, sim) in enumerate(zip(with_aa_reference, basal_reference, enzyme_conc)):
			row = i // cols
			col = i % cols
			ax = plt.subplot(gs[row, col])

			ax.plot(times, sim, label='Simulation')
			ax.plot([start, end], [with_aa, with_aa], 'r--', linewidth=1, label='Expected, with AA')
			ax.plot([start, end], [basal, basal], 'g--', linewidth=1, label='Expected, basal')

			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.set_title(aa_ids[i], fontsize=8)

			if row == rows - 1:
				ax.set_xlabel('Time (min)')
			if col == 0:
				ax.set_ylabel('Concentration (uM)')

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
