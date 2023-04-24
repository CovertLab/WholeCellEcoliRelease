"""
Compares expression of proteins containing rare arginines.
"""

import pickle
import os
import time

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

report_proteome = False
report_genes = False

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# Check variants
		variants = self.ap.get_variants()
		if not (0 in variants and 6 in variants):
			print('Early return -- Variants 0 and 6 of the'
				' trna_synthetase_kinetics variant is required for this plot.')
			return

		################################################################
		# Locate arginine codons across proteome
		arginine_codons = sim_data.relation.amino_acid_to_codons['ARG[c]']

		protein_to_composition = {}
		has_CGG_CGG = []
		has_AGA_AGA = []
		has_AGG_AGG = []
		has_no_arginine = []
		for protein, codon_sequence in sim_data.relation._codon_sequences.items():

			composition = '_'.join(codon_sequence)
			protein_to_composition[protein] = composition

			if 'CGG_CGG' in composition:
				has_CGG_CGG.append(protein)

			if 'AGA_AGA' in composition:
				has_AGA_AGA.append(protein)

			if 'AGG_AGG' in composition:
				has_AGG_AGG.append(protein)

			if np.all([codon not in composition for codon in arginine_codons]):
				has_no_arginine.append(protein)

		################################################################
		# Group proteins by if they have tandem CGGs, AGAs or AGGs
		set_CGG_CGG = set(has_CGG_CGG)
		set_AGA_AGA = set(has_AGA_AGA)
		set_AGG_AGG = set(has_AGG_AGG)
		all_proteins = set(sim_data.relation._codon_sequences.keys())
		set_no_arginine = set(has_no_arginine)
		other_proteins = (all_proteins
			- set_CGG_CGG - set_AGA_AGA - set_AGG_AGG - set_no_arginine)

		if report_proteome:
			print('Reporting number of proteins with the following rare'
				' arginine tandems:')
			print('-'*80)
			print(f'CGG_CGG:\t{len(set_CGG_CGG)}')
			print(f'AGA_AGA:\t{len(set_AGA_AGA)}')
			print(f'AGG_AGG:\t{len(set_AGG_AGG)}')
			set_CGG_AGA_intersection = set_CGG_CGG.intersection(set_AGA_AGA)

			print('Intersections')
			print(f'\tCGG, CGG and AGA, AGA:\t{len(set_CGG_AGA_intersection)}')
			print(f'\tAGA, AGA and AGG, AGG:\t'
				f'{len(set_AGA_AGA.intersection(set_AGG_AGG))}')
			print(f'\tAGG, AGG and CGG, CGG:\t'
				f'{len(set_AGG_AGG.intersection(set_CGG_CGG))}')
			print(f'\tAll 3:\t'
				f'{len(set_CGG_CGG.intersection(set_AGA_AGA, set_AGG_AGG))}')
			print(f'Proteins without arginine:\t{len(set_no_arginine)}')
			print(f'All other proteins:\t{len(other_proteins)}')

			argA = 'N-ACETYLTRANSFER-MONOMER'
			print(f'ArgA in CGG_CGG:\t{argA in set_CGG_CGG}')
			print(f'ArgA in AGA_AGA:\t{argA in set_AGA_AGA}')
			print(f'ArgA in AGG_AGG:\t{argA in set_AGG_AGG}')

		################################################################
		# Plot
		fig, ax = plt.subplots(1, 1, figsize=(8, 2.5))

		# Get expression
		variants = [0, 6]
		variant_to_expression = {}
		for variant in variants:
			cell_paths = self.ap.get_cells(
				variant=[variant], generation=np.arange(8, 12))
			if len(cell_paths) == 0:
				print('Early return -- this plots focuses on data generated'
					' between generations 8 and 11.')
				return
			expression = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts').mean(axis=0)
			variant_to_expression[variant] = expression

		# Get data
		data = []
		monomers = sim_data.process.translation.monomer_data['id'].tolist()
		for i, proteins in enumerate([set_CGG_CGG, set_AGA_AGA, set_AGG_AGG,
				other_proteins, set_no_arginine]):

			# Get protein compartments
			proteins = list(proteins)
			proteins_compartments = sim_data.getter.get_compartments(proteins)
			proteins_with_compartments = [
				f'{protein}[{loc[0]}]' for (protein, loc)
				in zip(proteins, proteins_compartments)]

			# Get fold changes
			fold_changes = []
			for protein in proteins_with_compartments:
				protein_index = monomers.index(protein)
				fold_changes.append(
					(1 + variant_to_expression[6][protein_index])
					/ (1 + variant_to_expression[0][protein_index]))

			# Record fold changes
			data.append(fold_changes)

			# Record CGG, CGG data
			if i == 0:
				CGG_CGG_proteins = proteins
				CGG_CGG_fold_changes = fold_changes

		# Report genes in CGG, CGG group
		if report_genes:
			for i in np.argsort(CGG_CGG_fold_changes)[::-1]:
				print(f'{CGG_CGG_proteins[i]}\t{CGG_CGG_fold_changes[i]}')

		# Plot
		color = '#6b6ecf'
		xticks = np.arange(5)
		violins = ax.violinplot(data, positions=xticks)
		for pc in violins['bodies']:
			pc.set_facecolor(color)
		for key in ['cmaxes', 'cmins', 'cbars']:
			violins[key].set_color(color)
			violins[key].set_linewidth(1.0)

		for x, fold_changes in enumerate(data):
			median = np.median(fold_changes)
			ax.plot([x, x], [median, median], marker='s', markersize=4,
				color=color, mew=0, linewidth=1.0)
			ax.text(x + 0.05, median, f'{median:.2f}',
				ha='left', va='center', fontsize=7)

		ax.set_xticks(xticks)
		ax.set_xticklabels(['CGG, CGG', 'AGA, AGA', 'AGG, AGG',
			'Other Arginine Codons', 'No Arginines'])
		ax.tick_params(axis='both', which='major', labelsize=7)

		ax.set_title('Expression of Proteins with Rare Arginine Codons',
			fontsize=11)
		ax.set_xlabel('Protein Features', fontsize=9)
		ax.set_ylabel('Expression Fold Change', fontsize=9)
		ax.set_yscale('log')

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		################################################################
		# Report statistics
		def calculate_Z_score(x, mu, sigma):
			Z = (x - mu) / sigma
			return Z

		# Reference distribution: "No Arginines"
		print('Reference distribution: No Arginines')
		reference = data[4]
		mu = np.median(reference)
		sigma = np.std(reference)

		# Tandem CGGs
		x = np.median(data[0])
		Z = calculate_Z_score(x, mu, sigma)
		p = 2 * scipy.stats.norm.sf(abs(Z))
		print(f'\tTandem CGGs:\tmedian = {x:.2f}, Z = {Z:.2f}, p-value = {p:.1e}')

		# Tandem AGAs
		x = np.median(data[1])
		Z = calculate_Z_score(x, mu, sigma)
		p = 2 * scipy.stats.norm.sf(abs(Z))
		print(f'\tTandem AGAs:\tmedian = {x:.2f}, Z = {Z:.2f}, p-value = {p:.1e}')

		# Tandem AGGs
		x = np.median(data[2])
		Z = calculate_Z_score(x, mu, sigma)
		p = 2 * scipy.stats.norm.sf(abs(Z))
		print(f'\tTandem AGGs:\tmedian = {x:.2f}, Z = {Z:.2f}, p-value = {p:.1e}')

		# Other arginine codons
		x = np.median(data[3])
		Z = calculate_Z_score(x, mu, sigma)
		p = 2 * scipy.stats.norm.sf(abs(Z))
		print(f'\tOther arginine codons:\tmedian = {x:.2f}, Z = {Z:.2f}, p-value = {p:.1e}')


if __name__ == "__main__":
	Plot().cli()
