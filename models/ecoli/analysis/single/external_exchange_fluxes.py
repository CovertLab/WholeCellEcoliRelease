from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile) as f:
			sim_data = cPickle.load(f)
		carbon_sources = sim_data.moleculeGroups.carbon_sources

		# Exchange flux
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		initial_time = main_reader.readAttribute('initialTime')
		time = (main_reader.readColumn('time')[1:] - initial_time) / 3600

		fba_results = TableReader(os.path.join(simOutDir, 'FBAResults'))
		ex_flux = fba_results.readColumn('externalExchangeFluxes')[1:, :]
		ex_molecules = fba_results.readAttribute('externalMoleculeIDs')
		mean_fluxes = np.mean(ex_flux, axis=0)

		# Find the top imports and exports to label
		n_top_exchanges = 5
		sorted_idx = np.argsort(mean_fluxes)
		highest = {ex_molecules[i] for i in sorted_idx[:n_top_exchanges]}
		highest.update({ex_molecules[i] for i in sorted_idx[-n_top_exchanges:]})

		# Create plot
		plt.figure(figsize=(8.5, 11))

		## Plot fluxes for specific molecules of interest
		molecule_ids = ['OXYGEN-MOLECULE[p]', 'CARBON-DIOXIDE[p]']
		n_molecules = len(molecule_ids)
		rows = n_molecules + 3  # 2 for export/import and 1 for carbon sources
		cols = 1

		# Plot carbon sources
		ax = plt.subplot(rows, cols, 1)
		for carbon_source in carbon_sources:
			if carbon_source not in ex_molecules:
				continue

			molecule_flux = -ex_flux[:, ex_molecules.index(carbon_source)]
			if np.any(molecule_flux > 0):
				ax.plot(time, molecule_flux, label=carbon_source)
		ylim = ax.get_ylim()
		ax.set_ylim(max(0, ylim[0]), ylim[1])
		ax.set_ylabel('Carbon source import\n(mmol/gDCW/hr)')
		ax.tick_params(labelsize=8, which='both', direction='out')
		ax.legend()

		# Plot other uptake of interest
		for index, molecule in enumerate(molecule_ids):
			if molecule not in ex_molecules:
				continue

			molecule_flux = ex_flux[:, ex_molecules.index(molecule)]
			ave_flux = np.average(molecule_flux)
			if ave_flux < 0:
				molecule_flux = -1 * molecule_flux
				ave_flux *= -1
				direction = 'Import'
			else:
				direction = 'Export'

			ax = plt.subplot(rows, cols, index + 2)
			ax.plot(time, molecule_flux)

			y_range = np.min([
				np.abs(np.max(molecule_flux) - ave_flux),
				np.abs(np.min(molecule_flux) - ave_flux)
				])
			ymin = np.floor(ave_flux - y_range)
			ymax = np.ceil(ave_flux + y_range)
			ax.set_ylim([ymin, ymax])

			abs_max = np.max(molecule_flux)
			abs_min = np.min(molecule_flux)

			ax.text(0, ymin + 0.1*(ymax - ymin), 'Max: {:.1f}\nMin: {:.1f}'
				.format(abs_max, abs_min), fontsize=8)

			ax.set_ylabel('{} {}\n(mmol/gDCW/hr)'.format(molecule, direction), fontsize=8)
			ax.set_title(molecule, fontsize=8)
			ax.tick_params(labelsize=8, which='both', direction='out')

		## Plot all fluxes
		ax_export = plt.subplot(rows, cols, n_molecules + 2)
		ax_import = plt.subplot(rows, cols, n_molecules + 3)
		for mol, mean, flux in sorted(zip(ex_molecules, mean_fluxes, ex_flux.T),
				key=lambda v: v[1]):
			# Adjust for directionality
			if mean > 0:
				ax = ax_export
			elif mean < 0:
				ax = ax_import
				flux = -1 * flux
			else:
				continue

			# Create label if one of highest fluxes
			if mol in highest:
				label = '{}: {:.1f}'.format(mol, flux.mean())
			else:
				label = None

			ax.plot(time, flux, label=label)

		## Plot formatting
		ax_export.set_yscale('symlog')
		ax_import.set_yscale('symlog')
		handles, labels = ax_export.get_legend_handles_labels()
		ax_export.legend(handles[::-1], labels[::-1], fontsize=6)
		ax_import.legend(fontsize=6)
		ax_export.set_title('Export Fluxes', fontsize=8)
		ax_import.set_title('Import Fluxes', fontsize=8)
		ax_export.tick_params(labelsize=8, which='both', direction='out')
		ax_import.tick_params(labelsize=8, which='both', direction='out')
		ax_import.set_xlabel('Time (hr)', fontsize=8)
		ax_export.set_ylabel('Export Flux\n(mmol/gDCW/hr)', fontsize=8)
		ax_import.set_ylabel('Import Flux\n(mmol/gDCW/hr)', fontsize=8)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
