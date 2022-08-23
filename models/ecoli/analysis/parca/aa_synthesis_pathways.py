"""
Shows a comparison of kcat and KI values associated with amino acid synthesis
pathways.  Useful as validation for the kcat data and identifying any issues
with KI values that might not be inconsistent with expected amino acid
concentrations.
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		aa_synthesis_pathways = metabolism.aa_synthesis_pathways
		conc = metabolism.concentration_updates.concentrations_based_on_nutrients

		# Load kcat data for allosteric reactions
		amino_acids = sorted(aa_synthesis_pathways)
		kcat_data = np.array([
			aa_synthesis_pathways[aa]['kcat_data'].asNumber(1 / units.s)
			for aa in amino_acids
			])
		kcat_calc = np.array([
			aa_synthesis_pathways[aa]['kcat'].asNumber(1 / units.s)
			for aa in amino_acids
			])

		# Get kcat limits for plot
		positive_mask = (kcat_data > 0) & (kcat_calc > 0)
		kcat_min = min(kcat_data[positive_mask].min(), kcat_calc[positive_mask].min())
		kcat_max = max(kcat_data.max(), kcat_calc.max())
		kcat_range = [kcat_min, kcat_max]

		# Statistics
		r, p = pearsonr(np.log10(kcat_data[positive_mask]), np.log10(kcat_calc[positive_mask]))
		n = np.sum(positive_mask)

		# Calculate expected inhibition in minimal media condition
		aa_inhibition = np.array([
			[0, 0] if aa_synthesis_pathways[aa]['ki'] is None else
			[1 / (1 + conc('minimal')[aa] / ki).asNumber() for ki in aa_synthesis_pathways[aa]['ki']]
			for aa in amino_acids
			])

		# Ranges for plot
		x_amino_acids = np.arange(len(amino_acids))

		# Plot data
		plt.figure(figsize=(4, 12))
		gs = gridspec.GridSpec(3, 1)

		## kcat scatter comparison
		plt.subplot(gs[0, 0])
		plt.loglog(kcat_data[positive_mask], kcat_calc[positive_mask], 'o', markeredgewidth=0, markersize=10, alpha=0.5)
		plt.loglog(kcat_range, kcat_range, 'k--', linewidth=1, alpha=0.5)
		plt.xlabel('kcat, from data')
		plt.ylabel('kcat, calculated')
		plt.title(f'r={r:.2f} (p={p:.2g}, n={n})', fontsize=8)
		self.remove_border()

		## kcat bar comparison
		plt.subplot(gs[1, 0])
		plt.bar(x_amino_acids, np.log2(kcat_calc / kcat_data))
		plt.axhline(0, linestyle='--', color='k')
		plt.xticks(x_amino_acids, amino_acids, rotation=45, fontsize=6, ha='right')
		plt.ylabel('log2(kcat, calculated / kcat, from data)')
		self.remove_border()

		## Inhibition comparison
		plt.subplot(gs[2, 0])
		plt.bar(x_amino_acids, aa_inhibition[:, 1], label='upper limit KI')
		plt.bar(x_amino_acids, aa_inhibition[:, 0], label='lower limit KI')
		plt.xticks(x_amino_acids, amino_acids, rotation=45, fontsize=6, ha='right')
		plt.ylabel('Fraction of max flux at AA conc')
		self.remove_border()
		plt.legend()

		## Save plot
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
