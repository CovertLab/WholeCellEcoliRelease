"""
Plot ribosome profile on argA mRNAs.
"""

import pickle
import os
import copy

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import amino_acid_to_label
from wholecell.utils import units


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		variants = self.ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		# Analyze variants 0 (control) and 6 (ArgRS)
		if not (0 in variants and 6 in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 6.')
			return

		# Codon sequence
		amino_acid = 'ARG[c]'
		monomer = 'N-ACETYLTRANSFER-MONOMER'
		codon_sequence = np.array(sim_data.relation._codon_sequences[monomer])

		# Identify codons of interest (coi)
		codons = sim_data.relation.amino_acid_to_codons[amino_acid]
		coi_indexes = np.zeros(codon_sequence.shape[0], dtype=np.bool)
		for codon in codons:
			coi_indexes[codon_sequence == codon] = True
		non_coi_indexes = np.logical_not(coi_indexes)

		# tRNAs
		rna_data = sim_data.process.transcription.rna_data
		free_trnas = rna_data['id'][rna_data['is_tRNA']]
		charged_trnas = sim_data.process.transcription.charged_trna_names

		# Mapping from codons to trnas
		codons_to_trnas = sim_data.relation.trnas_to_codons.T.astype(np.bool)
		all_codons = sim_data.relation.codons

		# Plot
		fig, axes = plt.subplots(3, 2,
			figsize=(3.75, 3),
			gridspec_kw={'height_ratios': [1, 2, 1]},
			)
		color_non_coi = '#e7e7e7'# '#d9d9d9'
		colors = [
			'#6b6ecf', # blue
			'#b5cf6b', # green
			'#e7ba52', # yellow
			'#ce6dbd', # pink
			]

		# Order of codons across top
		codons = ['CGU', 'CGG', 'CGC', 'AGG']

		# Order of codons plotted (layers)
		plot_ordered_codons = ['CGU', 'CGG', 'CGC', 'AGG']
		codon_to_colors = dict(zip(codons, colors))

		# Sequence guides
		guide_locations = np.linspace(
			np.where(codon_sequence == 'CGU')[0][0],
			np.where(codon_sequence == 'CGU')[0][-1],
			len(codons),
			)

		for ax in axes[0, :]:
			for codon in plot_ordered_codons:
				i = codons.index(codon)
				color = codon_to_colors[codon]
				guide_location = guide_locations[i]
				for x in np.where(codon_sequence == codon)[0]:
					ax.plot([x, x], [0, 1], linewidth=0.5, color=color)
					ax.plot([x, guide_location], [1, 2], linewidth=0.5, color=color)
				ax.plot([guide_location, guide_location], [2, 3], linewidth=0.5, color=color)
				ax.scatter(guide_location, 3, marker='o', color=color)

		for ax in axes[2, :]:
			guide_location = guide_locations.min() + (guide_locations.max() - guide_locations.min()) / 2
			for x in np.where(non_coi_indexes)[0]:
				ax.plot([x, x], [0, -1], linewidth=0.5, color=color_non_coi)
				ax.plot([x, guide_location], [-1, -2], linewidth=0.5, color=color_non_coi)
			ax.plot([guide_location, guide_location], [-2, -3], linewidth=0.5, color=color_non_coi)
			ax.scatter(guide_location, -3, marker='o', color=color_non_coi)

		# Plot data
		n_assume_feasible = 49
		for i, variant in enumerate([0, 6]):

			# Get cells
			cell_paths = self.ap.get_cells(variant=[variant])

			# Instantiate ribosomes data
			ribosomes = np.zeros(len(codon_sequence))

			# Consider time steps within each cell cycle
			for sim_dir in cell_paths:
				simOutDir = os.path.join(sim_dir, 'simOut')
				main_reader = TableReader(os.path.join(simOutDir, 'Main'))
				time_step = main_reader.readColumn('timeStepSec')

				trna_reader = TableReader(
					os.path.join(simOutDir, 'TrnaCharging'))
				ribosome_positions = trna_reader.readColumn(
					'ribosome_positions_argA')
				collision_removal = trna_reader.readColumn(
					'collision_removed_ribosomes_argA')
				ribosome_initiations = trna_reader.readColumn(
					'ribosome_initiation_argA')

				for j in range(2, ribosome_positions.shape[0]):

					# Get current ribosome positions
					current = []
					for k, n in enumerate(ribosome_positions[j]):
						for l in range(n):
							current.append(k)
					original_current = copy.deepcopy(current)

					# No active ribosomes
					if len(current) == 0:
						continue

					# Get current time step
					time = time_step[j]

					# Calculated number of newly initiated ribosomes
					n_initiations = np.logical_not(
						np.isnan(ribosome_initiations[j - 1])).sum()

					# Get previous ribosome positions
					previous = []
					for k, n in enumerate(ribosome_positions[j - 1]):
						for l in range(n):
							previous.append(k)

					previous_collision_removed = []
					for k, n in enumerate(collision_removal[j - 1]):
						for l in range(n):
							previous_collision_removed.append(k)
					for x in previous_collision_removed:
						if x in previous:
							k = previous.index(x)
							_ = previous.pop(k)
						elif x == 0:
							# This ribosome removed by a collision event
							# was a newly initiated ribosome -- so let
							# it be removed from the count of n_initiations
							n_initiations -= 1
						else:
							print('edge case 1: collision removed a non-existent ribosome?')

					# If an mRNA is recorded as having a new ribosome
					# initated on it, but the ribosome positions suggest
					# that the new ribosome did not processes (take any
					# steps), then consider the initiation was not
					# successful.
					if n_initiations == 1:
						if original_current == previous:
							n_initiations -= 1
						elif len(original_current) == len(previous) and 443 not in original_current:
							n_initiations -= 1

					# current = np.array(current)
					# previous = np.array(previous)

					# Newly initiated ribosomes start from the beginning
					for k in range(n_initiations):
						_current = current[0]
						steps = _current + 1
						t = time / steps
						ribosomes[0:_current + 1] += t
						current = current[1:]

					# For each ribosome, match its previous position to
					# its current position
					for _previous in previous[::-1]:

						# Skip terminating ribosomes
						if _previous == ribosome_positions.shape[1] - 1:
							continue

						# Assess potential current positions
						for _current in current[::-1]:
							steps = _current - _previous

							# Paused ribosomes
							if steps == 0:
								ribosomes[_previous] += time
								k = current.index(_current)
								_ = current.pop(k)
								break

							# Processed ribosomes
							elif steps <= n_assume_feasible:
								t = time / steps
								ribosomes[_previous + 1:_current + 1] += t
								k = current.index(_current)
								_ = current.pop(k)
								break

							# Other
							else:
								print('edge case 2: increase n_assume_feasible?')
								pass

			# Calculate frequency of observations
			ribosome_profile = ribosomes / ribosomes.sum()

			ax = axes[1, i]
			for x, y in zip(np.where(non_coi_indexes)[0], ribosome_profile[non_coi_indexes]):
				ax.plot([x, x], [0, y], color=color_non_coi, linewidth=1)

			for x, y in zip(np.where(coi_indexes)[0], ribosome_profile[coi_indexes]):
				color = codon_to_colors[codon_sequence[x]]
				ax.plot([x, x], [0, y], color=color, linewidth=1)

			ax = axes[0, i]
			total = ribosomes.sum()
			for j, codon in enumerate(codons):
				guide_location = guide_locations[j]
				indexes = np.where(codon_sequence == codon)[0]
				f = ribosomes[indexes].sum() / total
				ax.text(guide_location, 3.5, f'{f * 100:.2f}%\n{codon}', ha='center', va='bottom', fontsize=7)

			ax = axes[2, i]
			f = ribosomes[non_coi_indexes].sum() / total
			guide_location = guide_locations.min() + (guide_locations.max() - guide_locations.min()) / 2
			ax.text(guide_location, -3.5, f'Non-{amino_acid_to_label[amino_acid]} Codons\n{f * 100:.2f}%', ha='center', va='top', fontsize=7)

		# Format
		axes[1, 0].set_ylabel('Fraction of Time', fontsize=7)
		axes[1, 0].set_xlabel('Codon Sequence', fontsize=7)
		axes[1, 1].set_xlabel('Codon Sequence', fontsize=7)

		for ax in axes[0, :]:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['left'].set_visible(False)
			ax.set_xticks([])
			ax.set_yticks([])
			ax.set_ylim([0, 3.5])

		for ax in axes[1, :]:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.set_yscale('log')
			ax.set_yticks([1e-5, 1e-3, 1e-1])
			ax.set_xticks([0, 200, 442])
			ax.tick_params(axis='both', which='major', labelsize=7)

		for ax in axes[2, :]:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['left'].set_visible(False)
			ax.set_xticks([])
			ax.set_yticks([])
			ax.set_ylim([-3.5, 0.5])

		xlim = axes[1, 1].get_xlim()
		for ax in axes.reshape(-1):
			ax.set_xlim(xlim)

		axes[1, 0].set_ylim(axes[1, 1].get_ylim())

		plt.tight_layout(h_pad=0, w_pad=0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
