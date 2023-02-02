"""
Plots fold change in molecule concentrations.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
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
		atp = 'ATP[c]'
		amino_acid = 'ARG[c]'
		trnas = sim_data.relation.amino_acid_to_trnas[amino_acid]
		synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]

		# Get data
		def get_concentration(volume, n_molecule):
			c_molecule = (1
				/ sim_data.constants.n_avogadro
				/ volume
				* n_molecule
				).asNumber(units.umol / units.L)
			return c_molecule

		molecules = ['ArgRS', 'Arginine', 'Arginine tRNAs', 'ATP']
		molecules_label = ['ArgRS', 'Arginine', 'Arginine\ntRNAs', 'ATP']
		variant_0_molecule_to_data = {molecule: [] for molecule in molecules}
		variant_6_molecule_to_data = {molecule: [] for molecule in molecules}

		# Variant 0, control
		variant = 0
		seeds = self.ap.n_seed
		generations = self.ap.n_generation
		molecule_to_concentrations = {
			molecule: np.nan * np.ones((seeds, generations))
			for molecule in molecules}

		for seed in range(seeds):
			for generation in range(generations):
				sim_dir = self.ap.get_cells(
					variant=[variant], seed=[seed], generation=[generation])

				if len(sim_dir) == 0:
					break

				sim_dir = sim_dir[0]
				sim_out_dir = os.path.join(sim_dir, 'simOut')
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

				volume = (1e-15
					* units.L
					* mass_reader.readColumn('cellVolume')
					)
				(n_atp, n_amino_acid, n_trnas, n_synthetase,
					) = read_bulk_molecule_counts(
					sim_out_dir,
					([atp], [amino_acid], trnas, [synthetase],)
					)
				molecule_to_concentrations['ATP'][seed, generation]\
					= get_concentration(volume, n_atp).mean()
				molecule_to_concentrations['Arginine'][seed, generation]\
					= get_concentration(volume, n_amino_acid).mean()
				molecule_to_concentrations['Arginine tRNAs'][seed, generation]\
					= get_concentration(volume, n_trnas.sum(axis=1)).mean()
				molecule_to_concentrations['ArgRS'][seed, generation]\
					= get_concentration(volume, n_synthetase).mean()

		for seed in range(seeds):
			for generation in range(1, generations):
				for molecule in molecules:
					current = molecule_to_concentrations[molecule][seed, generation]

					if np.isnan(current):
						continue
					previous = molecule_to_concentrations[molecule][seed, generation - 1]
					variant_0_molecule_to_data[molecule].append(current / previous)


		# Variant 6, experiment
		variant = 6
		for seed in range(self.ap.n_seed):
			for generation in range(self.ap.n_generation):

				sim_dir = self.ap.get_cells(
					variant=[variant], seed=[seed], generation=[generation])
				if len(sim_dir) == 0:
					break
				sim_dir = sim_dir[0]

				sim_out_dir = os.path.join(sim_dir, 'simOut')
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
				protein_mass = mass_reader.readColumn('proteinMass')

				if protein_mass[-1] / protein_mass[0] < 1.1:

					# Current concentration
					volume = (1e-15
						* units.L
						* mass_reader.readColumn('cellVolume')
						)
					(n_atp, n_amino_acid, n_trnas, n_synthetase,
						) = read_bulk_molecule_counts(
						sim_out_dir,
						([atp], [amino_acid], trnas, [synthetase],)
						)
					c_atp_current = get_concentration(volume, n_atp)
					c_amino_acid_current = get_concentration(volume, n_amino_acid)
					c_trnas_current = get_concentration(volume, n_trnas.sum(axis=1))
					c_synthetase_current = get_concentration(volume, n_synthetase)

					# Previous concentration
					sim_dir = self.ap.get_cells(
						variant=[variant],
						seed=[seed],
						generation=[generation - 1]
						)[0]
					sim_out_dir = os.path.join(sim_dir, 'simOut')
					mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
					volume = (1e-15
						* units.L
						* mass_reader.readColumn('cellVolume')
						)
					(n_atp, n_amino_acid, n_trnas, n_synthetase,
						) = read_bulk_molecule_counts(
						sim_out_dir,
						([atp], [amino_acid], trnas, [synthetase],)
						)
					c_atp_previous = get_concentration(volume, n_atp)
					c_amino_acid_previous = get_concentration(volume, n_amino_acid)
					c_trnas_previous = get_concentration(volume, n_trnas.sum(axis=1))
					c_synthetase_previous = get_concentration(volume, n_synthetase)

					# Calculate fold change
					variant_6_molecule_to_data['ATP'].append(c_atp_current.mean() / c_atp_previous.mean())
					variant_6_molecule_to_data['Arginine'].append(c_amino_acid_current.mean() / c_amino_acid_previous.mean())
					variant_6_molecule_to_data['Arginine tRNAs'].append(c_trnas_current.mean() / c_trnas_previous.mean())
					variant_6_molecule_to_data['ArgRS'].append(c_synthetase_current.mean() / c_synthetase_previous.mean())

					# Break to next seed
					break

		# Print stats
		# z score: z = (x - mu) / sigma
		# where x = a raw score, and mu and sigma are the mean and
		# standard deviations of the population
		x = np.mean(variant_0_molecule_to_data['Arginine'])
		mu = np.mean(variant_6_molecule_to_data['Arginine'])
		sigma = np.std(variant_6_molecule_to_data['Arginine'])
		z = (x - mu) / sigma
		print(f'z score: {z:.1f}')
		print(f'control / experiment: {x/mu:.1f}')

		# Plot
		def format_violins(violins, color, color_accent):
			for pc in violins['bodies']:
				pc.set_facecolor(color)
				pc.set_edgecolor('none')
				pc.set_alpha(alpha)

			for key in ['cbars', 'cmaxes', 'cmins']:
				violins[key].set_color(color_accent)
				violins[key].set_linewidth(0.5)

		color_control = '#bdbdbd'
		color_experiment = '#6b6ecf'
		color_control_accent = '#636363'
		color_experiment_accent = '#393b79'
		alpha = 0.7
		width = 0.4

		fontsize_heading = 11
		fontsize_label = 9
		fontsize_legend = 7

		fig, ax = plt.subplots(1, 1, figsize=(2.5, 2))
		center_positions = np.arange(len(molecules))

		def get_bar_parameters(data):
			mean = np.nanmean(data, axis=1)
			std = np.nanstd(data, axis=1)
			return mean, std
		error_kw = {'lw': 0.5, 'capsize': 2, 'capthick': 0.5}

		# Variant 0, control
		data = [np.log2(variant_0_molecule_to_data[molecule]) for molecule in molecules]
		median, yerr = get_bar_parameters(data)
		ax.bar(center_positions, median,
			width=-width, align='edge', color=color_control,
			yerr=yerr, error_kw=error_kw)

		# Variant 6, experiment
		positions = center_positions + (width / 2)
		data = [np.log2(variant_6_molecule_to_data[molecule]) for molecule in molecules]
		median, yerr = get_bar_parameters(data)
		ax.bar(center_positions, median,
			width=width, align='edge', color=color_experiment,
			yerr=yerr, error_kw=error_kw)

		# Format
		ax.set_title('Similarity to Previous\nCell Cycles', fontsize=fontsize_heading)
		ax.set_xlabel('ArgRS and Substrates', fontsize=fontsize_label)
		ax.set_ylabel('log2(Fold Change\nBetween Cell Cycles)', fontsize=fontsize_legend)

		ax.set_xticks(np.arange(len(molecules)))
		ax.set_xticklabels(molecules_label, fontsize=fontsize_legend)
		ax.tick_params(axis='both', which='major', labelsize=fontsize_legend)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.margins(0.01)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
