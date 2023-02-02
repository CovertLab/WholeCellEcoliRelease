"""
Plots molecules and enzymes participating in the arginine biosynthesis
pathway

"""

import pickle
import os

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		variants = self.ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		# Analyze variants 0 (control) and 6 (ArgRS)
		if not (0 in variants and 6 in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 6.')
			return

		# Molecules
		molecules = [
			'GLT[c]',
			'L-CITRULLINE[c]',
			'ARG[c]',
			'N-ACETYLTRANSFER-CPLX[c]',
		]

		molecule_to_label = {
			'GLT[c]': 'Glutamate',
			'L-CITRULLINE[c]': 'Citrulline',
			'ARG[c]': 'Arginine',
			'N-ACETYLTRANSFER-CPLX[c]': 'N-Acetylglutamate Synthase',
		}

		reactions = [
			'N-ACETYLTRANSFER-RXN',
		]

		reaction_to_label = {
			'N-ACETYLTRANSFER-RXN': 'Reaction Flux',
			}

		plotting_order = [
			'GLT[c]',
			'N-ACETYLTRANSFER-RXN',
			'N-ACETYLTRANSFER-CPLX[c]',
			'L-CITRULLINE[c]',
			'ARG[c]',
			]

		# Get data
		variants = [0, 6]
		seeds = [8, 8]
		variant_to_time = {}
		variant_to_reactions = {}
		variant_to_molecules = {}
		division_events = []

		cell_density = sim_data.constants.cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
		mmol_per_g_per_h = units.mmol / units.g / units.h

		for i, (variant, seed) in enumerate(zip(variants, seeds)):
			cell_paths = self.ap.get_cells(variant=[variant], seed=[seed])

			# Get data
			time = read_stacked_columns(
				cell_paths, 'Main', 'time', remove_first=True).reshape(-1) / 60
			volume = 1e-15 * units.L * read_stacked_columns(
				cell_paths, 'Mass', 'cellVolume', remove_first=True)
			(n_molecules,) = read_stacked_bulk_molecules(
				cell_paths, (molecules,), remove_first=True)
			c_molecules = (1
				/ sim_data.constants.n_avogadro
				/ volume
				* n_molecules
				).asNumber(units.umol / units.L)
			variant_to_time[variant] = time
			variant_to_molecules[variant] = c_molecules

			cell_mass = read_stacked_columns(
				cell_paths, 'Mass', 'cellMass', remove_first=True).reshape(-1)
			dry_mass = read_stacked_columns(
				cell_paths, 'Mass', 'dryMass', remove_first=True).reshape(-1)
			raw_reaction_fluxes = read_stacked_columns(
				cell_paths, 'FBAResults', 'reactionFluxes', remove_first=True)
			coefficient = dry_mass / cell_mass * cell_density

			sim_out_dir = os.path.join(cell_paths[0], 'simOut')
			fba_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
			reaction_ids = fba_reader.readAttribute('reactionIDs')

			reaction_indexes = [reaction_ids.index(reaction) for reaction in reactions]
			reaction_fluxes = (
					COUNTS_UNITS
					/ MASS_UNITS
					/ TIME_UNITS
					* raw_reaction_fluxes[:, reaction_indexes]
					/ coefficient[:, None]).asNumber(mmol_per_g_per_h)
			variant_to_reactions[variant] = reaction_fluxes

			if variant == 6:
				for sim_dir in cell_paths:
					sim_out_dir = os.path.join(sim_dir, 'simOut')
					main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
					division_events.append(main_reader.readColumn('time')[0] / 60)

		# Calculate average arginine concentration
		arginine_index = molecules.index('ARG[c]')
		c_arg_control = variant_to_molecules[0][:, arginine_index].mean()

		# Identify period of time when arginine is depleted
		c_arg = variant_to_molecules[6][:, arginine_index]
		arg_depletion = c_arg < c_arg_control - 30 # uM
		i_start = np.where(arg_depletion)[0][0]
		i_end = np.where(arg_depletion)[0][-1]

		# Identify when the N-acetyl transfer reaction is 0
		reaction_index = reactions.index('N-ACETYLTRANSFER-RXN')
		reaction_flux = variant_to_reactions[6][:, reaction_index]
		is_zero = np.where(reaction_flux <= 0)[0]
		continuous_array = np.arange(is_zero[0], is_zero[-1] + 1)
		i_start_zero_rxn_flux = is_zero[0]
		i_end_zero_rxn_flux = is_zero[-1]

		# Plot
		fig, axes = plt.subplots(2, 5, figsize=(8, 1.75),
			gridspec_kw={'height_ratios': [1, 4]})

		color_control = '#bdbdbd'
		color_experiment = '#6b6ecf'
		color_highlight_red = '#d6616b' # red
		color_highlight_green = '#b5cf6b' # green
		color_highlight_yellow = '#e7ba52' # yellow
		color_highlight_pink = '#ce6dbd' # pink

		linewidth = 1

		colors = [color_control, color_experiment]
		fontsize_label = 9
		fontsize_legend = 7

		def annotate_reaction_flux(ax):
			x1 = variant_to_time[6][i_start_zero_rxn_flux]
			x2 = variant_to_time[6][i_end_zero_rxn_flux]
			y = 1
			ax.plot([x1, x2], [y, y], color=color_highlight_red)
			return

		def annotate_arginine_depletion(ax):
			x1 = variant_to_time[6][i_start]
			x2 = variant_to_time[6][i_end]
			y = 0
			ax.plot([x1, x2], [y, y], color=color_highlight_yellow)
			return

		################################################################
		window_time = 60 * 10 # seconds
		window_size = window_time // 2 # 2 second time steps
		window = np.ones(window_size)

		plot_zero_reaction_flux = False

		for i, molecule_or_reaction in enumerate(plotting_order):
			ax = axes[1, i]

			# Molecule concentration
			if molecule_or_reaction in molecules:
				molecule = molecule_or_reaction
				j = molecules.index(molecule)
				ylabel = 'uM'
				control_concentration = np.median(variant_to_molecules[0][:, j])
				metabolite_concentration = variant_to_molecules[6][:, j]

				if molecule in ['GLT[c]', 'ARG[c]', 'L-CITRULLINE[c]']:
					ylabel = 'mM'
					control_concentration /= 1000
					metabolite_concentration /= 1000

				ax.set_ylabel(ylabel, fontsize=fontsize_legend)
				ax.axhline(
					control_concentration,
					color=color_control,
					linewidth=linewidth,
					)

				smoothened = np.convolve(
					metabolite_concentration, window, mode='same') / window_size
				ax.plot(
					variant_to_time[6][window_size:-window_size],
					smoothened[window_size:-window_size],
					color=color_experiment,
					linewidth=linewidth,
					)

				# Annotate
				ax = axes[0, i]
				ax.set_title(molecule_to_label[molecule], fontsize=fontsize_label)
				if plot_zero_reaction_flux:
					annotate_reaction_flux(ax)
				annotate_arginine_depletion(ax)

			# Reaction flux
			elif molecule_or_reaction in reactions:
				reaction = molecule_or_reaction
				j = reactions.index(reaction)
				ax.axhline(
					np.median(variant_to_reactions[0][:, j]),
					color=color_control,
					linewidth=linewidth,
					)

				reaction_flux = variant_to_reactions[6][:, j]
				smoothened = np.convolve(
					reaction_flux, window, mode='same') / window_size
				ax.plot(
					variant_to_time[6][window_size:-window_size],
					smoothened[window_size:-window_size],
					color=color_experiment,
					linewidth=linewidth,
					)

				# Format
				ax.set_ylabel('mmol/gDCW/h', fontsize=fontsize_legend)
				ax.set_yscale('log')
				ax.set_yticks([1e-2, 1e-1, 1e0])
				ax.set_ylim([1e-2, 1e0])

				# Annotate
				ax = axes[0, i]
				ax.set_title(reaction_to_label[reaction], fontsize=fontsize_label)
				annotate_reaction_flux(ax)
				plot_zero_reaction_flux = True
				annotate_arginine_depletion(ax)

		################################################################

		# Format all axes
		generations = [0, 9, 14, 19]

		for ax in axes[0, :]:
			ax.set_ylim([-1, 2])
			ax.set_xlim(axes[1, 0].get_xlim())
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['left'].set_visible(False)
			ax.set_xticks([])
			ax.set_yticks([])

		for ax in axes[1, :]:
			ax.tick_params(axis='both', which='major', labelsize=fontsize_legend)
			ax.set_xticks(np.array(division_events)[generations])
			ax.set_xticklabels(generations)
			ax.set_xlabel('Generation', fontsize=fontsize_legend)
			ax.yaxis.offsetText.set_fontsize(fontsize_legend)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)


		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
