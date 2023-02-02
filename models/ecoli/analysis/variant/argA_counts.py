"""
Plots argA molecules.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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

		# Check for the required sim variants
		variant_control = 0
		variant_experiment = 6
		variants = [variant_control, variant_experiment]
		sim_variants = self.ap.get_variants()

		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		if False in [variant in sim_variants for variant in variants]:
			print('This analysis script is designed for the variant:'
				f' trna_synthetase_kinetics variant indexes: {variants}.')
			return

		# Molecules
		argA_monomer = 'N-ACETYLTRANSFER-MONOMER[c]'
		argA_complex = 'N-ACETYLTRANSFER-CPLX[c]'
		complexation_reaction = 'N-ACETYLTRANSFER-CPLX_RXN'

		# Load data and plot
		def plot(cell_paths, axes, color):
			for generation, sim_dir in enumerate(cell_paths):

				# Listeners used
				sim_out_dir = os.path.join(sim_dir, 'simOut')
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				ribosome_reader = TableReader(os.path.join(sim_out_dir, 'RibosomeData'))
				complexation_reader = TableReader(os.path.join(sim_out_dir, 'ComplexationListener'))
				monomer_deg_reader = TableReader(os.path.join(sim_out_dir, 'MonomerDegradationListener'))
				replication_reader = TableReader(os.path.join(sim_out_dir, 'ReplicationData'))

				# Load data
				time = main_reader.readColumn('time') / 60
				(n_complex, n_monomer,) = read_bulk_molecule_counts(
					sim_out_dir, ([argA_complex], [argA_monomer],))

				i = complexation_reader.readAttribute('reactionIDs').index(complexation_reaction)
				complexation_events = complexation_reader.readColumn('complexationEvents')[:, i]

				i = ribosome_reader.readAttribute('monomerIds').index(argA_monomer)
				monomer_synthesis_events = ribosome_reader.readColumn('n_terminated')[:, i]
				monomer_degradation_events = monomer_deg_reader.readColumn('monomers_degraded')[:, i]

				# Plot
				ax = axes[0]
				ax.plot(time, monomer_synthesis_events, color=color)

				ax = axes[1]
				ax.plot(time, monomer_degradation_events, color=color)

				ax = axes[2]
				ax.plot(time, n_monomer, color=color)

				ax = axes[3]
				ax.plot(time, complexation_events, color=color)

				ax = axes[4]
				ax.plot(time, n_complex, color=color)

				if generation == 0:
					print(n_complex[0])

				# Indicate 0 complexes
				if 0 in n_complex:
					if np.all(n_complex == 0):
						ax.plot([time[0], time[-1]], [90, 90], color_highlight_red)

				# Cell division events
				for ax in axes:
					ax.axvline(time[-1], color='k', linestyle='--',
						linewidth=0.2)

				# Text labels
				x = time[0] + (time[-1] - time[0]) / 2
				if generation % 2 == 0:
					if generation == 0:
						text = f'Gen\n0'
					else:
						text = f'{generation}'
					axes[0].text(x, 3.5, text, fontsize=fontsize_annotation,
						va='bottom', ha='center')
			return

		# Initialize plot
		fig, axes = plt.subplots(5, 2, figsize=(3, 3), sharey='row', sharex='col',
			gridspec_kw={'width_ratios': [1, 10]})

		color_highlight_red = '#d6616b' # red
		fontsize_annotation = 7
		fontsize_label = 7
		fontsize_heading = 9
		linewidth_dashed = 1

		# Plot control variant
		color_control = '#bdbdbd'
		cell_paths = self.ap.get_cells(variant=[0], seed=[8], generation=[0])
		plot(cell_paths, axes[:, 0], color_control)

		# Plot experimental variant
		color_experiment = '#6b6ecf'
		cell_paths = self.ap.get_cells(variant=[6], seed=[8])
		plot(cell_paths, axes[:, 1], color_experiment)

		# Format
		axes[0, 1].set_title('Monomer Synthesis', fontsize=fontsize_heading)
		axes[1, 1].set_title('Monomer Degradation', fontsize=fontsize_heading)
		axes[2, 1].set_title('Monomer', fontsize=fontsize_heading)
		axes[3, 1].set_title('Complexation', fontsize=fontsize_heading)
		axes[4, 1].set_title('Hexamer', fontsize=fontsize_heading)

		axes[0, 0].set_ylabel('Number\nof Events', fontsize=fontsize_label)
		axes[1, 0].set_ylabel('Number\nof Events', fontsize=fontsize_label)
		axes[2, 0].set_ylabel('Number of\nMolecules', fontsize=fontsize_label)
		axes[3, 0].set_ylabel('Number\nof Events', fontsize=fontsize_label)
		axes[4, 0].set_ylabel('Number of\nMolecules', fontsize=fontsize_label)
		axes[4, 0].set_xlabel('Time (min)', fontsize=fontsize_label)
		axes[4, 1].set_xlabel('Time (min)', fontsize=fontsize_label)

		# Set margins
		for ax in axes[:, 1]:
			ax.margins(0.01)

		# Monomer synthesis
		for ax in axes[0, :]:
			ax.set_ylim([0, 4.5])
			ax.set_yticks([0, 3])

		# Monomer degradation
		for ax in axes[1, :]:
			ax.set_yticks([0, 1])

		# Monomer counts
		for ax in axes[2, :]:
			ax.set_yticks([0, 6])

		# Complexation
		for ax in axes[3, :]:
			ax.set_yticks([0, 1])

		# Hexamer counts
		for ax in axes[4, :]:
			ax.set_yticks([0, 100])

		# Set ticks
		y_portion = 0.08
		for row in range(axes.shape[0]):
			ax = axes[row, 0]
			y_max = ax.get_ylim()[1]
			y_min = -(y_max * y_portion / (1 - y_portion))
			ax.set_ylim([y_min, y_max])

		for ax in axes.reshape(-1):
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(axis='both', which='major', labelsize=fontsize_annotation)

		plt.tight_layout(w_pad=0, h_pad=0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
