"""
Plots average count of trna synthetases.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		# Measurement counts
		schmidt_ids = validation_data.protein.schmidt2015Data['monomerId']
		schmidt_counts = validation_data.protein.schmidt2015Data['glucoseCounts']

		# Simulated counts
		model_ids = sim_data.process.translation.monomer_data['id']
		model_counts = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts')

		# Unify
		all_model_counts, all_schmidt_counts = get_simulated_validation_counts(
			schmidt_counts, model_counts, schmidt_ids, model_ids)

		# trna synthetases
		synthetases = []
		for amino_acid in sim_data.molecule_groups.amino_acids:
			if amino_acid == 'L-SELENOCYSTEINE[c]':
				continue
			synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]
			synthetases += sim_data.process.complexation.get_monomers(synthetase)['subunitIds'].tolist()

		model_ids_list = model_ids.tolist()
		indexes = [model_ids_list.index(molecule) for molecule in synthetases]
		highlight_model_counts, highlight_schmidt_counts = get_simulated_validation_counts(
			schmidt_counts, model_counts[:, indexes], schmidt_ids, model_ids[indexes])

		# Plot
		color = '#bdbdbd'
		color_highlight = '#6b6ecf'

		fig, ax = plt.subplots(1, 1, figsize=(2.7, 2.7))

		# All proteins
		x = np.log10(all_schmidt_counts + 1)
		y = np.log10(all_model_counts + 1)
		label = 'All Proteins'

		# Counts >= 30
		threshold = np.log10(30)
		mask = np.logical_and(x >= threshold, y >= threshold)
		r, _ = pearsonr(x[mask], y[mask])
		r_squared = r**2
		label += f'\nCopies >= 30: R2 = {r_squared:.2f} (n = {sum(mask)})'

		# Counts < 30
		mask = np.logical_not(mask)
		r, _ = pearsonr(x[mask], y[mask])
		r_squared = r**2
		label += f'\nCopies < 30: R2 = {r_squared:.2f} (n = {sum(mask)})'

		ax.scatter(x, y, s=8, marker='.', linewidths=0, color=color, label=label)

		# Highlight trna synthetases
		x = np.log10(highlight_schmidt_counts + 1)
		y = np.log10(highlight_model_counts + 1)
		r, _ = pearsonr(x, y)
		r_squared = r**2
		ax.scatter(x, y, s=8, marker='o', edgecolors=color_highlight,
			facecolors='none', label=f'tRNA Synthetases: R^2 = {r_squared:.2f} (n = {len(x)}))')

		# Format
		ax.set_xlabel('Measured\nlog10(Number of Molecules + 1)', fontsize=9)
		ax.set_ylabel('Simulated\nlog10(Number of Molecules + 1)', fontsize=9)
		ax.legend(loc='upper left', fontsize=7)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(axis='both', which='major', labelsize=7)

		# Guide lines
		ax_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
		ax_max = np.ceil(max(ax.get_xlim()[1], ax.get_ylim()[1]))
		ax.plot([ax_min, ax_max], [ax_min, ax_max], color='k', linewidth=1)
		ax.plot([ax_min, ax_max - 1], [ax_min + 1, ax_max], color='k', linewidth=1, linestyle='--')
		ax.plot([ax_min + 1, ax_max], [ax_min, ax_max - 1], color='k', linewidth=1, linestyle='--')
		ax.set_xlim([ax_min, ax_max])
		ax.set_ylim([ax_min, ax_max])
		ax.set_xticks([0, 1, 2, 3, 4, 5, 6])
		ax.set_yticks([0, 1, 2, 3, 4, 5, 6])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
