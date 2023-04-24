"""
Locates position of ribosomes along the argA mRNA when they experienced
premature termination.
"""

import pickle
import os

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		# Locate premature terminations
		raw_premature_terminations = read_stacked_columns(
			cell_paths, 'TrnaCharging', 'collision_removed_ribosomes_argA')
		premature_terminations = raw_premature_terminations.sum(axis=0)

		# Plot
		fig, axes = plt.subplots(2, 1, figsize=(8, 3.2))

		# Premature ribosome termination
		ax = axes[0]
		color = '#636363'
		percent = premature_terminations / premature_terminations.sum() * 100
		ax.plot(np.arange(percent.shape[0]), percent, color=color)

		# Highlight the 2 peaks
		color_highlight = '#b5cf6b'
		x0 = 22
		x1 = 25
		ax.plot(np.arange(x0, x1), percent[x0:x1], color=color_highlight)

		x0 = 152
		x1 = 155
		ax.plot(np.arange(x0, x1), percent[x0:x1], color=color_highlight)

		# Label
		y = percent[23]
		ax.text(23, y, f'{y: .1f}%', va='bottom', ha='center', fontsize=7)

		y = percent[153] + percent[154]
		ax.text(153.5, y, f'{y: .1f}%', va='bottom', ha='center', fontsize=7)

		# Indicate CGGs
		y = 110
		dy = 8
		ax.scatter([23, 153.5], y * np.ones(2), marker='v', linewidth=0, color=color_highlight)
		ax.text(23, y + dy, 'CGG', va='bottom', ha='center', fontsize=7)
		ax.text(153.5, y + dy, 'CGG, CGG', va='bottom', ha='center', fontsize=7)

		# Format
		ax.set_xticks([0, 100, 200, 300, 400, 442])
		ax.set_title('Premature Ribosome Terminations', fontsize=11)
		ax.set_xlabel('argA Codon Sequence', fontsize=9)
		ax.set_ylabel('Percent of Terminations', fontsize=9)
		ax.set_ylim([ax.get_ylim()[0], 130])
		ax.set_xlim([-10, 452])
		ax.tick_params(axis='both', which='major', labelsize=7)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		################################################################
		# All arginine codons
		ax = axes[1]

		argA = 'N-ACETYLTRANSFER-MONOMER'
		codon_sequence = sim_data.relation._codon_sequences[argA]
		arginine_codons = sim_data.relation.amino_acid_to_codons['ARG[c]']

		codon_to_color = {
			'CGU': '#6b6ecf', # blue
			'CGG': '#b5cf6b', # green
			'CGC': '#e7ba52', # yellow
			'AGG': '#ce6dbd', # pink
			}

		y = 1.3
		dy = 0.16
		for x, codon in enumerate(codon_sequence):
			if codon in arginine_codons:
				color = codon_to_color[codon]
				ax.plot([x, x], [0, 1], color=color)

				# Highlight tandem arginine codons
				previous_codon = codon_sequence[x - 1]
				if previous_codon in arginine_codons:
					previous_color = codon_to_color[previous_codon]
					ax.scatter(x - 1, y, marker='v', linewidth=1,
						color='none', edgecolors=previous_color)
					ax.scatter(x, y, marker='v', linewidth=1,
						color='none', edgecolors=color)

					if x > 280:
						text = f'{previous_codon},\n{codon}'
					else:
						text = f'{previous_codon}, {codon}'
					ax.text(x - 0.5, y + dy, text,
						va='bottom', ha='center', fontsize=7)

		# Legend
		lines = []
		labels = []
		for codon, color in codon_to_color.items():
			lines.append(Line2D([0], [0], color=color))
			labels.append(codon)
		ax.legend(lines, labels, fontsize=7)

		ax.set_ylim([-0.08, 1.8])
		ax.set_xlim(axes[0].get_xlim())
		ax.set_xticks([0, 100, 200, 300, 400, 442])
		ax.set_yticks([])
		ax.tick_params(axis='both', which='major', labelsize=7)
		ax.set_xlabel('argA Codon Sequence', fontsize=9)
		ax.set_ylabel('Arginine Codons', fontsize=9)
		ax.set_title('Arginine Codon Locations', fontsize=11)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
