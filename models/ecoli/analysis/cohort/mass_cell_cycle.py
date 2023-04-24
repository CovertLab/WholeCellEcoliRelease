"""
Plots accumulation of mass over the cell cycle.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

print_stats = True

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells(generation=[0, 1])

		# Identify a representative cell closest to the average doubling time
		doubling_times = []
		for sim_dir in cell_paths:
			sim_out_dir = os.path.join(sim_dir, 'simOut')
			main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
			time = main_reader.readColumn('time') / 60
			doubling_times.append(time[-1] - time[0])
		doubling_times = np.array(doubling_times)
		mean_doubling_time = doubling_times.mean()
		i_representative = np.argsort(np.abs(doubling_times - mean_doubling_time))[0]
		
		sim_dir = cell_paths[i_representative]
		seed = int(sim_dir.split('/')[-3])
		generation = int(sim_dir.split('/')[-2].split('_')[1])
		print(f'Representative: seed {seed}, generation {generation}')

		# Retrive mass data of representative cell
		sim_dir = cell_paths[i_representative]
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		main_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
		time = main_reader.readColumn('time') / 60
		representative_time = time - time[0]
		mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
		representative_protein = mass_reader.readColumn('proteinMass')
		representative_dry = mass_reader.readColumn('dryMass')
		representative_non_protein = np.vstack((
				mass_reader.readColumn('dnaMass'),
				mass_reader.readColumn('rnaMass'),
				mass_reader.readColumn('smallMoleculeMass'),
				)).sum(axis=0)

		# Retrive mass fold change
		protein_mass_fold_change = []
		non_protein_mass_fold_change = []
		dry_mass_fold_change = []

		for sim_dir in cell_paths:
			sim_out_dir = os.path.join(sim_dir, 'simOut')
			mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

			# Protein
			protein = mass_reader.readColumn('proteinMass')
			protein_mass_fold_change.append(protein[-1] / protein[0])

			# Dry
			dry = mass_reader.readColumn('dryMass')
			dry_mass_fold_change.append(dry[-1] / dry[0])

			# Non-protein
			non_protein = np.vstack((
				mass_reader.readColumn('dnaMass'),
				mass_reader.readColumn('rnaMass'),
				mass_reader.readColumn('smallMoleculeMass'),
				)).sum(axis=0)
			non_protein_mass_fold_change.append(non_protein[-1] / non_protein[0])

		# z-score: (mean - population_mean) / std
		population_mean = 2 # expectation of doubling mass

		print('z scores')
		samples = [protein_mass_fold_change, non_protein_mass_fold_change, dry_mass_fold_change]
		names = ['Protein', 'Non-protein', 'Total (dry)']
		for sample, name in zip(samples, names):
			mean = np.mean(sample)
			std = np.std(sample)

			z = (mean - population_mean) / std
			print(f'\t{name}:\t{z:.2f}')

		# Plot
		def format_violins(violins, color, color_accent):
			for pc in violins['bodies']:
				pc.set_facecolor(color)
				pc.set_edgecolor('none')
				pc.set_alpha(alpha)

			for key in ['cbars', 'cmaxes', 'cmins']:
				violins[key].set_color(color_accent)
				violins[key].set_linewidth(0.5)

		fig, axes = plt.subplots(1, 2, figsize=(5, 2.5))
		guide_width = 0.4
		color_total = '#bdbdbd'
		color_protein = '#6b6ecf'
		color_non_protein = '#e7ba52'
		color_total_accent = '#636363'
		color_protein_accent = '#393b79'
		color_non_protein_accent = '#8c6d31'
		alpha = 0.7

		ax = axes[0]
		ax.plot(representative_time,
			representative_dry / representative_dry[0],
			color=color_total, label='Total (Dry)')
		ax.plot(representative_time,
			representative_protein / representative_protein[0],
			color=color_protein, label='Protein')
		ax.plot(representative_time,
			representative_non_protein / representative_non_protein[0],
			color=color_non_protein, label='Non-Protein')
		ax.set_xlabel('Time (min)', fontsize=9)
		ax.set_ylabel('Mass (Relative to Intial Mass)', fontsize=9)
		ax.set_title('Representative Cell Growth', fontsize=9)
		ax.legend(loc='best', fontsize=7)
		ax.set_yticks([1, 2, 3])

		ax = axes[1]
		violins = ax.violinplot(dry_mass_fold_change, positions=[0])
		format_violins(violins, color_total, color_total_accent)

		violins = ax.violinplot(protein_mass_fold_change, positions=[1])
		format_violins(violins, color_protein, color_protein_accent)

		violins = ax.violinplot(non_protein_mass_fold_change, positions=[2])
		format_violins(violins, color_non_protein, color_non_protein_accent)

		ax.set_yticks([1.5, 2, 2.5, 3.0])
		ax.set_xticks([0, 1, 2])
		ax.set_xticklabels(['Total (Dry)', 'Protein', 'Non-Protein'])
		ax.set_xlabel('Cellular Components', fontsize=9)
		ax.set_ylabel('Mass (Relative to Intial Mass)', fontsize=9)
		ax.set_title('Mass Fold Change', fontsize=9)

		# Format
		for ax in axes:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.axhline(2, color='k', linestyle='--', linewidth=1)
			ax.tick_params(axis='both', which='major', labelsize=7)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		if print_stats:
			print('Representative')
			print(f'\tDry:\t{representative_dry[-1] / representative_dry[0]:.3f}')
			print(f'\tProtein:\t{representative_protein[-1] / representative_protein[0]:.3f}')
			print(f'\tNon-Protein:\t{representative_non_protein[-1] / representative_non_protein[0]:.3f}')


			print('Summary of all cells')
			print(f'\tDry:\t{np.mean(dry_mass_fold_change):.3f}')
			print(f'\tProtein:\t{np.mean(protein_mass_fold_change):.3f}')
			print(f'\tNon-Protein:\t{np.mean(non_protein_mass_fold_change):.3f}')


if __name__ == '__main__':
	Plot().cli()
