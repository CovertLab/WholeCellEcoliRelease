"""
Analysis for metabolism_secretion_penalty variant to show impact of secretion
penalty on growth, flux correlation and exchanges.

TODO: read and plot central carbon validation correlation
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.metabolism_secretion_penalty import SECRETION_PENALTY
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


def remove_border(ax):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.set_xticks([])
	ax.tick_params(axis='y', labelsize=6)

def show_xaxis(ax, x_vals):
	ax.spines['bottom'].set_visible(True)
	ax.set_xticks(x_vals)
	ax.set_xticklabels(np.array(SECRETION_PENALTY)[x_vals], rotation=45, fontsize=6)

def plot_fluxes(gs, col, selected, x_vals, fluxes, labels):
	for i, idx in enumerate(selected):
		ax = plt.subplot(gs[i, col])
		data = fluxes[:, idx]
		if data.mean() < 0:
			label = 'Import'
			direction = -1
		else:
			label = 'Export'
			direction = 1
		ax.bar(x_vals, direction * data)
		ax.set_title(labels[idx][:30], fontsize=8)
		ax.set_ylabel('{} Flux\n(mmol/g DCW/hr)'.format(label), fontsize=8)
		remove_border(ax)
	show_xaxis(ax, x_vals)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		all_doubling_times = np.zeros(n_variants)
		all_exchange_fluxes = None
		ex_molecules = []
		x_vals = []

		for i, variant in enumerate(variants):
			if variant > len(SECRETION_PENALTY):
				continue
			x_vals.append(variant)

			doubling_times = []
			exchange_fluxes = []
			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				main_reader = TableReader(os.path.join(simOutDir, 'Main'))
				fba_reader = TableReader(os.path.join(simOutDir, 'FBAResults'))

				# Load data
				time = main_reader.readColumn('time')
				ex_molecules = fba_reader.readAttribute('externalMoleculeIDs')
				ex_flux = fba_reader.readColumn('externalExchangeFluxes')[1:, :]

				# Save cell data for this variant
				doubling_times.append((time[-1] - time[0]) / 60)
				exchange_fluxes.append(ex_flux)

			# Aggregate all data for this variant
			all_doubling_times[i] = np.mean(doubling_times)
			if all_exchange_fluxes is None:
				all_exchange_fluxes = np.zeros((n_variants, len(ex_molecules)))
			all_exchange_fluxes[i, :] = np.mean(np.vstack(exchange_fluxes), axis=0)

		# Determine number of plots to show
		n_data = 1  # data axes to plot below (not including fluxes)
		n_fluxes = 5  # number of high and low fluxes to plot
		sort_idx = np.argsort(all_exchange_fluxes.mean(axis=0))
		lowest_idx = sort_idx[:n_fluxes]
		highest_idx = sort_idx[::-1][:n_fluxes]

		# Create bar plots
		plt.figure(figsize=(8.5, 11))
		gs = gridspec.GridSpec(max(n_fluxes, n_data), 3)

		## Doubling time
		ax = plt.subplot(gs[0, 0])
		ax.bar(x_vals, all_doubling_times)
		ax.set_ylabel('Doubling Time\n(min)', fontsize=8)
		remove_border(ax)
		show_xaxis(ax, x_vals)

		## Lowest exchange fluxes
		col = 1
		plot_fluxes(gs, col, lowest_idx, x_vals, all_exchange_fluxes, ex_molecules)

		## Highest exchange fluxes
		col = 2
		plot_fluxes(gs, col, highest_idx, x_vals, all_exchange_fluxes, ex_molecules)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
