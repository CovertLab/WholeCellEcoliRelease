"""
Comparison of calculated amino acid uptake rates to measured uptake rates
from Zampieri et al. 2019.  Useful as validation of calculated values.
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.dataclasses.process.metabolism import DRY_MASS_UNITS, K_CAT_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


FLUX_UNITS = units.mmol / units.g / units.h


def subplot(x, y, x_err, labels, y_label, title):
	# Statistics
	r, p = pearsonr(np.log10(x), np.log10(y))
	n = len(labels)

	## Plot data
	min_rate = min(x.min(), y.min())
	max_rate = max(x.max(), y.max())
	plt.errorbar(x, y, xerr=x_err, fmt='o')
	plt.plot([min_rate, max_rate], [min_rate, max_rate], '--k')

	## Show point labels
	for aa, x, y in zip(labels, x, y):
		plt.text(x, 1.1 * y, aa, ha='center', fontsize=6)

	## Format axes
	plt.xlabel('Zampieri et al. max uptake flux\n(mmol/g DCW/hr)')
	plt.ylabel(y_label)
	plt.xscale('log')
	plt.yscale('log')
	plt.title(f'{title}\nr={r:.2f} (p={p:.2g}, n={n})', fontsize=8)

class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validation_data_file, 'rb') as f:
			validation_data = pickle.load(f)

		# Attributes from sim_data
		wcm_aa_ids = sim_data.molecule_groups.amino_acids
		specific_import_rates = sim_data.process.metabolism.specific_import_rates
		aa_supply_rates = sim_data.translation_supply_rate['minimal_plus_amino_acids']

		# Load validation data
		aa_ids = []
		val_rates = []
		val_lb = []
		val_ub = []
		for aa, rates in validation_data.amino_acid_uptake_rates.items():
			aa_ids.append(aa)
			val_rates.append(rates['uptake'].asNumber(FLUX_UNITS))
			val_lb.append(rates['LB'].asNumber(FLUX_UNITS))
			val_ub.append(rates['UB'].asNumber(FLUX_UNITS))
		val_rates = np.array(val_rates)
		val_error = np.vstack((val_rates - np.array(val_lb), np.array(val_ub) - val_rates))

		# Load WCM calculated data
		wcm_rates = []
		wcm_supply = []
		wcm_data = {
			aa[:-3]: rates
			for aa, rates in zip(wcm_aa_ids, specific_import_rates)
			}
		supply_data = {
			aa[:-3]: rates
			for aa, rates in zip(wcm_aa_ids, aa_supply_rates)
			}
		wcm_units = units.count / DRY_MASS_UNITS * K_CAT_UNITS  # TODO: store this in sim_data?
		for aa in aa_ids:
			wcm_rates.append((wcm_units * wcm_data[aa]).asNumber(FLUX_UNITS))
			wcm_supply.append(supply_data[aa].asNumber(FLUX_UNITS))
		wcm_rates = np.array(wcm_rates)
		wcm_supply = np.array(wcm_supply)

		# Create plot
		plt.figure(figsize=(5, 10))

		plt.subplot(2, 1, 1)
		subplot(val_rates, wcm_rates, val_error, aa_ids,
			'WCM uptake flux\n(mmol/g DCW/hr)',
			'Amino acid import rate validation')

		plt.subplot(2, 1, 2)
		subplot(val_rates, wcm_supply, val_error, aa_ids,
			'WCM translation supply in rich media\n(mmol/g DCW/hr)',
			'Translation vs import comparison')

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
