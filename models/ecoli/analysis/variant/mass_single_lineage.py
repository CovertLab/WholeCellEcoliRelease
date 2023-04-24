"""
Plots mass across cell components, and compares with a control.
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


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		variants = self.ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		variant_experiment = 6

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		if not (0 in variants and variant_experiment in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 6.')
			return

		def plot(cell_paths, axes, color):
			color_highlight = '#b5cf6b'
			found_dna = False
			found_rna = False
			found_protein = False
			found_total_dry = False

			for i, sim_dir in enumerate(cell_paths):
				# Listeners used
				sim_out_dir = os.path.join(sim_dir, 'simOut')
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

				# Load data
				time = main_reader.readColumn('time') / 60
				dna_mass = mass_reader.readColumn('dnaMass')
				rna_mass = mass_reader.readColumn('rnaMass')
				protein_mass = mass_reader.readColumn('proteinMass')
				dry_mass = mass_reader.readColumn('dryMass')

				# Plot
				ax = axes[0]
				if not found_dna and dna_mass[-1] / dna_mass[0] < 1.1:
					ax.plot(time, dna_mass, color=color_highlight, linewidth=linewidth)
					found_dna = True
				else:
					ax.plot(time, dna_mass, color=color, linewidth=linewidth)

				ax = axes[1]
				if not found_rna and rna_mass[-1] / rna_mass[0] < 1.1:
					ax.plot(time, rna_mass, color=color_highlight, linewidth=linewidth)
					found_rna = True
				else:
					ax.plot(time, rna_mass, color=color, linewidth=linewidth)

				ax = axes[2]
				if not found_protein and protein_mass[-1] / protein_mass[0] < 1.1:
					ax.plot(time, protein_mass, color=color_highlight, linewidth=linewidth)
					found_protein = True
				else:
					ax.plot(time, protein_mass, color=color, linewidth=linewidth)

				ax = axes[3]
				if not found_total_dry and dry_mass[-1] / dry_mass[0] < 1.1:
					ax.plot(time, dry_mass, color=color_highlight, linewidth=linewidth)
					found_total_dry = True
				else:
					ax.plot(time, dry_mass, color=color, linewidth=linewidth)

				# Cell division events
				for ax in axes:
					ax.axvline(
						time[-1], color='k', linestyle='--', linewidth=0.2)

				# Text labels
				x = time[0] + (time[-1] - time[0]) / 2
				y = 15

				if i % 2 == 0:
					if i == 0:
						text = f'Gen\n0'
					else:
						text = f'{i}'
					axes[0].text(x, y, text, fontsize=7,
						va='bottom', ha='center')
			return

		# Initialize plot
		fig, axes = plt.subplots(
			4, 2,
			figsize=(5, 2.5),
			gridspec_kw={'width_ratios': [1, 20]},
			sharey='row',
			sharex='col',
			)
		linewidth = 1

		# Plot mass of control variant
		color_control = '#969696'
		cell_paths = self.ap.get_cells(variant=[0], seed=[8], generation=[0])
		plot(cell_paths, axes[:, 0], color_control)

		# Plot mass of experimental variant
		color_experiment = '#6b6ecf'
		cell_paths = self.ap.get_cells(variant=[variant_experiment], seed=[8])
		plot(cell_paths, axes[:, 1], color_experiment)

		# Format
		# axes[0, 1].set_title('Mass of Cellular Components', fontsize=9)
		axes[0, 1].set_title('DNA', fontsize=9)
		axes[1, 1].set_title('RNA', fontsize=9)
		axes[2, 1].set_title('Protein', fontsize=9)
		axes[3, 1].set_title('Total Dry', fontsize=9)

		for ax in axes[:, 0]:
			ax.set_ylabel('fg', fontsize=7)

		axes[3, 0].set_xlabel('Time (min)', fontsize=7)
		axes[3, 1].set_xlabel('Time (min)', fontsize=7)

		# Control, x ticks
		for ax in axes[:, 0]:
			ax.set_xticks([0, 50])

		# Format
		for ax in axes[:, 1]:
			ax.margins(0.01) # default 0.05

		y_portion = 0.08
		# DNA
		for ax in axes[0, :]:
			ax.set_yticks([0, 20])
			y_max = int(np.ceil(ax.get_ylim()[1]))
			y_min = -int(np.ceil(y_max * y_portion / (1 - y_portion)))
			ax.set_ylim([y_min, y_max])

		# RNA
		for ax in axes[1, :]:
			ax.set_yticks([0, 200])
			y_max = int(np.ceil(ax.get_ylim()[1]))
			y_min = -int(np.ceil(y_max * y_portion / (1 - y_portion)))
			ax.set_ylim([y_min, y_max])

		# Protein
		for ax in axes[2, :]:
			ax.set_yticks([0, 300])
			y_max = int(np.ceil(ax.get_ylim()[1]))
			y_min = -int(np.ceil(y_max * y_portion / (1 - y_portion)))
			ax.set_ylim([y_min, y_max])

		# Dry Mass
		for ax in axes[3, :]:
			ax.set_yticks([0, 700])
			y_max = int(np.ceil(ax.get_ylim()[1]))
			y_min = -int(np.ceil(y_max * y_portion / (1 - y_portion)))
			ax.set_ylim([y_min, y_max])

		for ax in axes.reshape(-1):
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(axis='both', which='major', labelsize=7)

		plt.tight_layout(w_pad=0, h_pad=0) # default 1.08
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
