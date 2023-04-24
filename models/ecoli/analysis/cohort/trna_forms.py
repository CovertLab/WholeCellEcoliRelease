"""
Plots trna forms.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.plotting_tools import trna_to_label, amino_acid_to_label


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Aminoacyl-trna abundance
		rna_data = sim_data.process.transcription.rna_data
		raw_trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		charged_trnas = sim_data.process.transcription.charged_trna_names

		# Remove selC-tRNA
		selC = raw_trnas.index('selC-tRNA[c]')
		_ = raw_trnas.pop(selC)
		_ = charged_trnas.pop(selC)

		# Order alphabetically for plotting
		trnas = [trna_to_label(trna) for trna in raw_trnas]
		indexes = np.argsort(trnas)
		trnas = np.array(trnas)[indexes]
		raw_trnas = np.array(raw_trnas)[indexes]
		charged_trnas = np.array(charged_trnas)[indexes]

		# Get number of molecules
		(n_charged_trnas, n_free_trnas) = read_stacked_bulk_molecules(
			cell_paths, (charged_trnas, raw_trnas), remove_first=True)

		# Add spaces between amino acid famililes
		xticklabels = []
		previous = 'ala'
		for trna in trnas:
			amino_acid = trna[:3]

			if amino_acid != previous:
				xticklabels.append('')
				previous = amino_acid

			xticklabels.append(trna)
		xticklabels.append('')
		positions = np.where(np.array(xticklabels) != '')[0]

		# Plot
		fig, axes = plt.subplots(2, 1, figsize=(10, 3),
			gridspec_kw={'height_ratios': [9, 1]})
		color_free = '#bdbdbd'
		color_charged = '#9c9ede'
		color_free_accent = '#636363'
		color_charged_accent = '#393b79'
		color_free_edge = '#969696'
		color_charged_edge = '#5254a3'

		alpha = 0.5
		widths = 0.8

		ax = axes[0]
		violins_free = ax.violinplot(n_free_trnas, positions=positions, widths=widths)
		violins_charged = ax.violinplot(n_charged_trnas, positions=positions, widths=widths)
		ax.set_ylabel('Number of Molecules / Cell', fontsize=9)
		ax.set_yscale('log')

		# Legend
		ax.legend(handles=[
			mpatches.Patch(facecolor=color_charged, edgecolor=color_charged_edge, linewidth=0.5, alpha=alpha, label='Aminoacylated'),
			mpatches.Patch(facecolor=color_free, edgecolor=color_free_edge, linewidth=0.5, alpha=alpha, label='Unaminoacylated'),
			],
			bbox_to_anchor=(0.5, 1.05), loc='center', ncol=2, fontsize=7)

		# Format
		def format_violins(violins, color, color_accent, color_edge, facecolor):
			for pc in violins['bodies']:
				if facecolor:
					pc.set_facecolor(color)
				else:
					pc.set_facecolor('none')
				pc.set_edgecolor(color_edge)
				pc.set_linewidths(0.5)
				pc.set_alpha(alpha)

			for key in ['cbars', 'cmaxes', 'cmins']:
				violins[key].set_color(color_accent)
				violins[key].set_linewidth(0.5)

		format_violins(violins_free, color_free, color_free_accent, color_free_edge, facecolor=False)
		format_violins(violins_charged, color_charged, color_charged_accent, color_charged_edge, facecolor=True)

		ax.set_xticks(positions)
		ax.set_xticklabels(trnas, rotation=90)
		ax.tick_params(axis='both', which='major', labelsize=7)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# Unify x axes
		for ax in axes:
			ax.set_xlim([-1, max(positions) + 1])

		# Delineate amino acid families
		ax = axes[1]
		for spine in ax.spines.values():
			spine.set_visible(False)
		ax.set_xticks([])
		ax.set_yticks([])

		space = 0.4
		for amino_acid in sim_data.molecule_groups.amino_acids:
			if amino_acid == 'L-SELENOCYSTEINE[c]':
				continue

			_trnas = [trna_to_label(trna)
				for trna in sim_data.relation.amino_acid_to_trnas[amino_acid]]
			x = [xticklabels.index(trna) for trna in _trnas]
			x_min = min(x) - space
			x_max = max(x) + space
			ax.plot([x_min, x_max], [1, 1], color='k', linewidth=1)

			if len(_trnas) <= 3 or amino_acid == 'GLN[c]':
				text = _trnas[0][0].upper() + _trnas[0][1:3]
			else:
				text = amino_acid_to_label[amino_acid]
			ax.text(np.mean(x), 0.95, text, fontsize=7, ha='center', va='top')

		# Print stats
		print('Average (across all isoacceptors, across both forms):')
		all_trnas = np.hstack((n_free_trnas[:, indexes], n_charged_trnas[:, indexes]))
		print(f'\t{all_trnas.mean():.1f} molecules / cell')

		print('Lowest average:')
		all_trnas_mean = all_trnas.mean(axis=0)
		i = np.argmin(all_trnas_mean)
		x = all_trnas_mean[i]
		if i < len(trnas):
			label = f'{trnas[i]} (free)'
		else:
			label = f'{trnas[i - len(trnas)]} (charged)'
		print(f'\t{label}: {x:.1f}')

		print('Highest average:')
		i = np.argmax(all_trnas_mean)
		x = all_trnas_mean[i]
		if i < len(trnas):
			label = f'{trnas[i]} (free)'
		else:
			label = f'{trnas[i - len(trnas)]} (charged)'
		print(f'\t{label}: {x:.1f}')

		print('Standard deviations')
		all_trnas_std = all_trnas.std(axis=0)
		i = np.argmin(all_trnas_std)
		x = all_trnas_std[i]
		if i < len(trnas):
			label = f'{trnas[i]} (free)'
		else:
			label = f'{trnas[i - len(trnas)]} (charged)'
		print(f'\tMin:\t{label}: {x:.3f}')

		i = np.argmax(all_trnas_std)
		x = all_trnas_std[i]
		if i < len(trnas):
			label = f'{trnas[i]} (free)'
		else:
			label = f'{trnas[i - len(trnas)]} (charged)'
		print(f'\tMax:\t{label}: {x:.3f}')

		print('Standard deviations / mean')
		all_trnas_std = all_trnas.std(axis=0)
		relative_std = all_trnas_std / all_trnas_mean
		i = np.argmin(relative_std)
		x = relative_std[i]
		if i < len(trnas):
			label = f'{trnas[i]} (free)'
		else:
			label = f'{trnas[i - len(trnas)]} (charged)'
		print(f'\tMin:\t{label}: {x:.3f}')

		i = np.argmax(relative_std)
		x = relative_std[i]
		if i < len(trnas):
			label = f'{trnas[i]} (free)'
		else:
			label = f'{trnas[i - len(trnas)]} (charged)'
		print(f'\tMax:\t{label}: {x:.3f}')

		print('Aminoacylated fraction')
		fraction = n_charged_trnas / (n_charged_trnas + n_free_trnas)
		fraction_average = fraction.mean(axis=0)

		i = np.argmin(fraction_average)
		x = fraction_average[i]
		label = f'{trnas[i]}'
		print(f'\tMin:\t{label}: {x:.4f}')

		i = np.argmax(fraction_average)
		x = fraction_average[i]
		label = f'{trnas[i]}'
		print(f'\tMax:\t{label}: {x:.4f}')

		all_trnas_ranges = all_trnas.max(axis=0) - all_trnas.min(axis=0)
		relative_range = all_trnas_ranges / all_trnas_mean

		plt.tight_layout(h_pad=0, w_pad=0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
