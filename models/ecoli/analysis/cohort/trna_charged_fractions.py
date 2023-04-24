"""
Plots trna charged fractions up until the generation that lacks protein
mass accumulation
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
# from wholecell.analysis.plotting_tools import amino_acid_to_label
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

amino_acid_to_label = {
	'L-ALPHA-ALANINE[c]': 'Ala',
	'ARG[c]': 'Arg',
	'ASN[c]': 'Asn',
	'L-ASPARTATE[c]': 'Asp',
	'CYS[c]': 'Cys',
	'GLT[c]': 'Glu', # Glutamic acid
	'GLN[c]': 'Gln', # Glutamine
	'GLY[c]': 'Gly',
	'HIS[c]': 'His',
	'ILE[c]': 'Ile',
	'LEU[c]': 'Leu',
	'LYS[c]': 'Lys',
	'MET[c]': 'Met',
	'PHE[c]': 'Phe',
	'PRO[c]': 'Pro',
	'SER[c]': 'Ser',
	'THR[c]': 'Thr',
	'TRP[c]': 'Trp',
	'TYR[c]': 'Tyr',
	'VAL[c]': 'Val',
}

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# For each lineage, identify when protein mass stops accumulating
		cell_paths_to_analyze = []
		for seed in range(self.ap.n_seed):
			for generation in range(self.ap.n_generation):
				sim_dir = self.ap.get_cells(seed=[seed], generation=[generation])

				if len(sim_dir) == 0:
					break

				sim_out_dir = os.path.join(sim_dir[0], 'simOut')
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
				protein_mass = mass_reader.readColumn('proteinMass')

				if protein_mass[-1] / protein_mass[0] < 1.1:
					break

				else:
					cell_paths_to_analyze.append(sim_dir[0])

		print(f'Cells analyzed: {len(cell_paths_to_analyze)}')

		# Molecules
		rna_data = sim_data.process.transcription.rna_data
		free_trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		charged_trnas = sim_data.process.transcription.charged_trna_names

		amino_acids = sim_data.molecule_groups.amino_acids
		i = amino_acids.index('L-SELENOCYSTEINE[c]')
		_ = amino_acids.pop(i)
		amino_acids = np.array(amino_acids)

		# Get trnas
		(n_free_trnas, n_charged_trnas,) = read_stacked_bulk_molecules(
			cell_paths_to_analyze, (free_trnas, charged_trnas,), remove_first=True)

		# Group trnas by amino acid
		data = []
		for amino_acid in amino_acids:
			trnas = sim_data.relation.amino_acid_to_trnas[amino_acid]
			trna_indexes = [free_trnas.index(trna) for trna in trnas]
			free_trnas_sum = n_free_trnas[:, trna_indexes].sum(axis=1)
			charged_trnas_sum = n_charged_trnas[:, trna_indexes].sum(axis=1)
			fraction = charged_trnas_sum / (charged_trnas_sum + free_trnas_sum)
			data.append(fraction)

		# Convert from fraction to percent
		data = np.array(data) * 100

		# Report arginine, lowest, and highest
		means = np.mean(data, axis=1)
		stds = np.std(data, axis=1)
		mask = amino_acids != 'ARG[c]'

		i = np.where(amino_acids == 'ARG[c]')[0][0]
		print(f'ARG:\t{means[i]:.1f}% charged (ARG[c])')
		print(f'\trange: {data[i, :].min():.1f} to {data[i, :].max():.1f}')

		i = np.argmin(means[mask])
		print(f'Low:\t{means[mask][i]:.1f}% charged ({amino_acids[mask][i]})')

		i = np.argmax(means[mask])
		print(f'High:\t{means[mask][i]:.1f}% charged ({amino_acids[mask][i]})')

		# Plot
		fig, ax = plt.subplots(1, 1, figsize=(2.75, 2))
		fontsize_heading = 11
		fontsize_label = 9
		fontsize_legend = 7
		color = '#6b6ecf'
		color_accent = '#393b79'

		for i in range(len(amino_acids)):
			violin = ax.violinplot(data[i, :], positions=[i], widths=1.0)
			for pc in violin['bodies']:
				pc.set_facecolor(color)
				pc.set_edgecolor('none')
				pc.set_alpha(0.7)

			for key in ['cbars', 'cmaxes', 'cmins']:
				violin[key].set_color(color_accent)
				violin[key].set_linewidth(0.5)

		ax.set_title('tRNA\nAminoacylation', fontsize=fontsize_heading)
		ax.set_xlabel('Amino Acid Family', fontsize=fontsize_label)
		ax.set_ylabel('Aminoacylated\nPercent', fontsize=fontsize_legend)

		ax.set_xticks(np.arange(len(amino_acids)))
		ax.set_xticklabels(
			[amino_acid_to_label[amino_acid] for amino_acid in amino_acids],
			fontsize=fontsize_legend, rotation=90)

		ax.set_yticks([0, 50, 100])
		ax.set_yticklabels(['0%', '50%', '100%'])
		ax.tick_params(axis='both', which='major', labelsize=fontsize_legend)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.margins(0.01)

		y_portion = 0.08
		y_max = ax.get_ylim()[1]
		y_min = -(y_max * y_portion / (1 - y_portion))
		ax.set_ylim([y_min, y_max])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
