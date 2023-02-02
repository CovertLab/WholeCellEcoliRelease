"""
Compares variants of trna synthetase kinetics variant.
"""

import pickle
import os
import re

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import pearsonr

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Bremer & Dennis. EcoSal Plus. 2008
BnD_x = [100, 60, 40, 30, 24]
BnD_y = np.array([74, 139, 241, 408, 571]) * 1e3

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		# Analyze variants 0 (control) and 1 (previous)
		if not (0 in variants and 1 in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 1.')
			return

		# Get data
		variants = [0, 1]
		variant_to_proteome = {}
		variant_to_fluxome = {}

		mmol_per_g_per_h = units.mmol / units.g / units.h
		cell_density = sim_data.constants.cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

		# Identify shared reactions between Toya and WCM
		sim_dir = ap.get_cells()[0]
		fba_reader = TableReader(os.path.join(sim_dir, 'simOut', 'FBAResults'))
		wcm_reactions = np.array(fba_reader.readAttribute('reactionIDs'))
		toya_reactions = validation_data.reactionFlux.toya2010fluxes['reactionID']

		reactions = []
		for reaction in toya_reactions:
			if np.any([reaction in X for X in wcm_reactions]):
				reactions.append(reaction)

		for variant in variants:
			cell_paths = ap.get_cells(variant=[variant])

			# Proteome
			n_monomers = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts')
			variant_to_proteome[variant] = n_monomers.mean(axis=0)

			# Fluxome
			cell_mass = read_stacked_columns(cell_paths, 'Mass', 'cellMass')
			dry_mass = read_stacked_columns(cell_paths, 'Mass', 'dryMass')
			coefficient = dry_mass / cell_mass * cell_density

			reaction_fluxes = (
				COUNTS_UNITS
				/ MASS_UNITS
				/ TIME_UNITS
				* read_stacked_columns(cell_paths, 'FBAResults', 'reactionFluxes')
				/ coefficient).asNumber(mmol_per_g_per_h)

			fluxes = []
			for reaction in reactions:
				indexes = np.where([reaction in X for X in wcm_reactions])[0]
				is_reverse = ['(reverse)' in X for X in wcm_reactions[indexes]]

				if np.any(is_reverse):
					reverse_indexes = indexes[is_reverse]
					not_reverse_indexes = indexes[np.logical_not(is_reverse)]
					flux = (reaction_fluxes[:, not_reverse_indexes].sum(axis=1)
						+ (-1 * reaction_fluxes[:, reverse_indexes].sum(axis=1)))
				else:
					flux = reaction_fluxes[:, indexes].sum(axis=1)
				fluxes.append(flux.mean())

			variant_to_fluxome[variant] = fluxes

		################################################################
		# Get cells
		cell_paths = ap.get_cells(variant=[0])

		# Get measurements
		measurements = validation_data.dong1996
		measurement_growth_rates = (measurements['growth_rates']
			).asNumber(1 / units.h)

		# Get data
		simulation_growth_rate = (1
			/ sim_data.condition_to_doubling_time[sim_data.condition]
			.asNumber(units.h))

		rna_data = sim_data.process.transcription.rna_data
		free_trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		charged_trnas = sim_data.process.transcription.charged_trna_names

		i = free_trnas.index('selC-tRNA[c]')
		_ = free_trnas.pop(i)
		_ = charged_trnas.pop(i)

		(n_free_trnas, n_charged_trnas) = read_stacked_bulk_molecules(
			cell_paths, (free_trnas, charged_trnas,))
		n_trnas = n_free_trnas + n_charged_trnas
		volume = 1e-15 * units.L * read_stacked_columns(
			cell_paths, 'Mass', 'cellVolume')
		c_trnas = (1
			/ sim_data.constants.n_avogadro
			/ volume
			* n_trnas
			).asNumber(units.umol / units.L)

		# Compare simulated and measured tRNA concentrations
		x_values = []
		y_values = []
		y_stds = []

		for i, trna in enumerate(free_trnas):

			# Measurement
			x_array = measurements[trna]
			f = interp1d(measurement_growth_rates, x_array, copy=False)
			x_values.append(f([simulation_growth_rate])[0])			

			# Simulation
			y = c_trnas[:, i]
			y_values.append(y.mean())
			y_stds.append(y.std())

		# Total tRNA mass
		rna_data = sim_data.process.transcription.rna_data
		all_free_trnas = rna_data['id'][rna_data['is_tRNA']]
		all_charged_trnas = sim_data.process.transcription.charged_trna_names
		all_trnas = all_free_trnas.tolist() + all_charged_trnas
		(n_trnas,) = read_stacked_bulk_molecules(cell_paths, (all_trnas,))
		n_trnas = n_trnas.sum(axis=1)

		# Doubling time
		doubling_times = []
		for sim_dir in cell_paths:
			sim_out_dir = os.path.join(sim_dir, 'simOut')
			main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
			time = main_reader.readColumn('time') / 60
			doubling_times.append(time[-1] - time[0])

		################################################################

		# Plot
		width = 8.2
		fig, axes = plt.subplots(1, 4, figsize=(width, width / 4))
		color_proteome = '#969696' #'#bdbdbd'
		color_fluxome = '#969696'
		index_control = 0
		index_previous = 1

		################################################################
		ax = axes[0]
		# ax.set_title('Proteome', fontsize=9)
		ax.set_xlabel('Previous Model'
			# '\n Number of Molecules + 1'
			, fontsize=9)
		ax.set_ylabel('Current Model'
			# '\n Number of Molecules + 1'
			, fontsize=9)
		x = 1 + variant_to_proteome[index_previous]
		y = 1 + variant_to_proteome[index_control]
		ax.scatter(x, y, color=color_proteome, s=2)

		r, _ = pearsonr(x, y)
		r_squared = r**2
		ax.text(0.5, 1, f'R2 = {r_squared:.4f} (n = {len(x)})',
			transform=ax.transAxes, ha='center', va='top', fontsize=7)

		ax.set_xscale('log')
		ax.set_yscale('log')

		################################################################
		ax = axes[1]
		# ax.set_title('Fluxome', fontsize=9)
		ax.set_xlabel('Previous Model'
			# '\nmmol / g / h'
			, fontsize=9)
		ax.set_ylabel('Current Model'
			# '\nmmol / g / h'
			, fontsize=9)
		x = variant_to_fluxome[index_previous]
		y = variant_to_fluxome[index_control]
		ax.scatter(x, y, color=color_fluxome, s=4)

		r, _ = pearsonr(x, y)
		r_squared = r**2
		ax.text(0.5, 1, f'R2 = {r_squared:.4f} (n = {len(x)})',
			transform=ax.transAxes, ha='center', va='top', fontsize=7)

		################################################################
		ax = axes[2]
		ax.scatter(
			x_values,
			y_values,
			s=4,
			color=color_fluxome,
			)

		# Format
		ax.set_xlabel('Measurement\nConcentration (uM)', fontsize=9)
		ax.set_ylabel('Current Model\nConcentration (uM)', fontsize=9)

		ax.tick_params(axis='both', which='major', labelsize=7)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		ax_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
		ax_max = max(ax.get_xlim()[1], ax.get_ylim()[1])

		r, _ = pearsonr(x_values, y_values)
		r_squared = r**2
		# ax.set_title('tRNA Abundance', fontsize=9)
		ax.text(0.5, 1, f'R2 = {r_squared:.4f} (n = {len(free_trnas)})',
			transform=ax.transAxes, ha='center', va='top', fontsize=7)

		################################################################
		ax = axes[3]

		# Dennis & Bremer, 2008
		ax.plot(BnD_x, BnD_y,
			# color='#bdbdbd',
			color='#969696',
			marker='o', markersize=4,
			linestyle='--', zorder=0, label='Bremer & Dennis, 2008')

		# Whole cell model
		violins = ax.violinplot(n_trnas, positions=[np.mean(doubling_times)],
			widths=6)
		for pc in violins['bodies']:
			pc.set_facecolor('#6b6ecf')
			pc.set_edgecolor('none')
			pc.set_alpha(0.7)

		for key in ['cbars', 'cmaxes', 'cmins']:
			violins[key].set_color('#393b79')
			violins[key].set_linewidth(0.5)

		# # Legend
		# handles, labels = ax.get_legend_handles_labels()
		# ax.legend(handles= handles
		# 	+ [mpatches.Patch(color='#393b79', alpha=0.7, label='Simulation')],
		# 	loc='upper right',
		# 	fontsize=7,
		# 	)

		# Format
		# ax.set_title('tRNA Abundance\n', fontsize=9)
		ax.set_xlabel('Doubling Time (min)', fontsize=9)
		ax.set_ylabel('Molecules per Cell', fontsize=9)
		ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
		ax.yaxis.offsetText.set_fontsize(7)

		################################################################
		for ax in axes:
			if ax in axes[:3]:
				ax_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
				ax_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
				ax.set_xlim([ax_min, ax_max])
				ax.set_ylim([ax_min, ax_max])
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(axis='both', which='major', labelsize=7)

		ax = axes[1]
		ax.set_xticks([-5, 0, 5, 10])
		ax.set_yticks([-5, 0, 5, 10])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
