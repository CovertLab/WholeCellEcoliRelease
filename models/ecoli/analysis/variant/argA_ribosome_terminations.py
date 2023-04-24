"""
Compares successful and unsuccessful translation terminations.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import ConnectionPatch
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

		# Analyze variants 0 (control) and 4 (ArgRS drop)
		variants = list(set(variants).intersection({0, 6}))

		# Arginine codon locations
		dnaN_monomer = 'N-ACETYLTRANSFER-MONOMER'
		codon_sequence = np.array(sim_data.relation._codon_sequences[dnaN_monomer])
		arg_codon_mask = np.zeros(codon_sequence.shape[0], dtype=np.bool)
		for codon in sim_data.relation.amino_acid_to_codons['ARG[c]']:
			arg_codon_mask[np.where(codon_sequence == codon)[0]] = True

		variant_to_data = {}
		for x, variant in enumerate(variants):
			cell_paths = self.ap.get_cells(variant=[variant])

			# Load data
			unsuccessful = read_stacked_columns(
				cell_paths, 'TrnaCharging', 'collision_removed_ribosomes_argA')
			n_unsuccessful = unsuccessful.sum()

			total = 0
			for sim_dir in cell_paths:
				sim_out_dir = os.path.join(sim_dir, 'simOut')
				reader = TableReader(os.path.join(sim_out_dir, 'TrnaCharging'))
				initiated = reader.readColumn('ribosome_initiation_argA')
				total += np.sum(np.logical_not(np.isnan(initiated)))

			n_successful = total - n_unsuccessful
			n_unsuccessful_arg = unsuccessful.sum(axis=0)[:-1][arg_codon_mask].sum()
			n_unsuccessful_non_arg = n_unsuccessful - n_unsuccessful_arg
			variant_to_data[variant] = {
				'n_successful': n_successful,
				'n_unsuccessful': n_unsuccessful,
				'n_unsuccessful_arg': n_unsuccessful_arg,
				'n_unsuccessful_non_arg': n_unsuccessful_non_arg,
			}

		# Plot
		fig, axes = plt.subplots(2, 1, figsize=(1.5, 3))
		color_completed = '#b5cf6b'# '#6b6ecf'
		color_interrupted_arg = '#e7ba52'
		color_interrupted_non_arg = '#bdbdbd'

		blue = '#6b6ecf'
		blue_light = '#9c9ede'
		blue_dark = '#5254a3'
		gray = '#e7e7e7'
		white = '#ffffff'

		# width = 0.6
		# tick_size = width / 10
		outer_labels = ['Interrupted', 'Completed']
		outer_colors = [blue, gray]
		inner_labels = ['Arginine', 'Non-Argnine', '']
		inner_colors = [blue_dark, blue_light, white]

		width = 0.3
		wedgeprops = {'width': width, 'edgecolor': white, 'linewidth': 0.5}

		def plot_pie(ax, variant):

			# Outer pie
			outer_ratios = np.array([
				variant_to_data[variant]['n_unsuccessful'],
				variant_to_data[variant]['n_successful'],
			])
			print(f'Completed:\t{outer_ratios[1] / outer_ratios.sum() * 100: .1f}%')
			print(f'Interrupted:\t{outer_ratios[0] / outer_ratios.sum() * 100: .1f}%')

			angle = - outer_ratios[0] / outer_ratios.sum() * 360 / 2

			_ = ax.pie(
				outer_ratios,
				radius=1,
				startangle=angle,
				colors=outer_colors,
				wedgeprops=wedgeprops,
				)

			# Inner pie
			inner_ratios = np.array([
				variant_to_data[variant]['n_unsuccessful_arg'],
				variant_to_data[variant]['n_unsuccessful_non_arg'],
				variant_to_data[variant]['n_successful'],
			])
			print(f'\tArg:\t{inner_ratios[0] / outer_ratios[0] * 100: .3f}%')
			print(f'\tNon-Arg:\t{inner_ratios[1] / outer_ratios[0] * 100: .3f}%')

			# # Less than 1% does not visually appear, so set to 1%
			if variant == 0 and inner_ratios[0] / inner_ratios.sum() < 1:
				a = 0.01 * inner_ratios.sum()
				b = inner_ratios[0] + inner_ratios[1] - a
				inner_ratios[0] = a
				inner_ratios[1] = b

			_ = ax.pie(
				inner_ratios,
				radius=1-width,
				startangle=angle,
				colors=inner_colors,
				wedgeprops=wedgeprops,
				)

			return

		print('Optimized')
		ax = axes[0]
		ax.set_title('Optimized ArgRS kcat', fontsize=9)
		plot_pie(ax, 0)

		ax = axes[1]
		ax.set_title('Measured ArgRS kcat', fontsize=9)
		plot_pie(ax, 6)

		plt.tight_layout(w_pad=0, h_pad=0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
