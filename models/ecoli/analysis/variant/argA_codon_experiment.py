"""
Compares argA expression under different experimental codon replacements.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import ConnectionPatch
# noinspection PyUnresolvedReferences
import scipy.stats
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

		# Get expression of argA
		variants = self.ap.get_variants()
		if not np.all([variant in variants for variant in [0, 8, 6, 7]]):
			print('Exiting -- requires four specific variants.')
			return

		variants = [0, 8, 6, 7]
		monomers = sim_data.process.translation.monomer_data['id'].tolist()
		argA_monomer = 'N-ACETYLTRANSFER-MONOMER[c]'
		argA_index = monomers.index(argA_monomer)

		# Plot
		fig, ax = plt.subplots(1, 1, figsize=(8, 2.5))
		color_blue = '#6b6ecf'
		color_gray = '#636363'

		n_cells = []
		means = []
		stds = []
		n_data_points = []

		for x, variant in enumerate(variants):
			if variant == 6:
				cell_paths = self.ap.get_cells(variant=[variant],
					generation=np.arange(8, 12))
			else:
				cell_paths = self.ap.get_cells(variant=[variant])

			if variant in [0, 8]:
				color = color_gray
			else:
				color = color_blue
			n_cells.append(len(cell_paths))
			expression = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts')[:, argA_index]

			parts = ax.violinplot(expression, positions=[x])
			mean = np.mean(expression)

			ax.plot([x, x], [mean, mean], marker='s', markersize=4,
				color=color, mew=0, linewidth=1.0)
			ax.text(x + 0.05, mean, f'{mean:.1f}',
				ha='left', va='center', fontsize=7)

			for pc in parts['bodies']:
				pc.set_facecolor(color)
			for key in ['cmaxes', 'cmins', 'cbars']:
				parts[key].set_color(color)
				parts[key].set_linewidth(1)

			means.append(mean)
			stds.append(np.std(expression))
			n_data_points.append(len(expression))

		xticklabels = [
			'CGG\nOptimized ArgRS',
			'CGU\nOptimized ArgRS',
			'CGG\nMeasured ArgRS',
			'CGU\nMeasured ArgRS',
			]
		for i, n in enumerate(n_cells):
			xticklabels[i] += f'\n(n = {n} cells)'

		ax.set_xticks(np.arange(len(xticklabels)))
		ax.set_xticklabels(xticklabels)
		ax.tick_params(axis='both', which='major', labelsize=7)
		ax.set_xlabel('Perturbation Experiments', fontsize=9)
		ax.set_ylabel('Number of Molecules per Cell', fontsize=9)
		ax.set_title('ArgA Expression', fontsize=11)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		################################################################
		# Report statistics
		def calculate_Z_score(variant_reference, variant_test):
			i_reference = variants.index(variant_reference)
			i_test = variants.index(variant_test)
			
			mu = means[i_reference]
			sigma = stds[i_reference]
			x = means[i_test]

			Z = (x - mu) / sigma
			return Z

		# Optimized ArgRS, CGG (v0) vs CGU (v8)
		Z = calculate_Z_score(0, 8)
		p = 2 * scipy.stats.norm.sf(abs(Z))

		print(f'Optimized ArgRS, CGG vs CGU: Z = {Z:.2f}, p-value = {p:.1e}')

		# Measured ArgRS, CGG (v6) vs CGU (v7)
		Z = calculate_Z_score(6, 7)
		p = 2 * scipy.stats.norm.sf(abs(Z))
		print(f'Measured ArgRS, CGG vs CGU: Z = {Z:.2f}, p-value = {p:.1e}')


if __name__ == "__main__":
	Plot().cli()
