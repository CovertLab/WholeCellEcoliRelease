"""
Plots impact of changing the argRS kcat.
"""

import pickle
import os
import copy

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


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
		amino_acid = 'ARG[c]'
		enzyme_monomer = 'N-ACETYLTRANSFER-MONOMER[c]'
		enzyme_complex = 'N-ACETYLTRANSFER-CPLX[c]'
		reaction = 'N-ACETYLTRANSFER-RXN'
		mrna = 'EG10063_RNA[c]'

		rna_data = sim_data.process.transcription.rna_data
		all_free_trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		all_charged_trnas = np.array(
			sim_data.process.transcription.charged_trna_names)

		free_trnas = sim_data.relation.amino_acid_to_trnas[amino_acid]
		trna_indexes = [all_free_trnas.index(trna) for trna in free_trnas]
		charged_trnas = all_charged_trnas[trna_indexes].tolist()

		synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]

		molecules = [amino_acid, synthetase, enzyme_monomer, enzyme_complex]\
			+ free_trnas + charged_trnas

		mrna_ids = rna_data['id'][rna_data['is_mRNA']].tolist()
		mrna_index = mrna_ids.index(mrna)

		codon_sequence = np.array(
			sim_data.relation._codon_sequences[enzyme_monomer.split('[c]')[0]])
		codons = sim_data.relation.amino_acid_to_codons[amino_acid]
		codon_to_position = {codon: [] for codon in codons}
		codon_to_position['other'] = []
		for i, codon in enumerate(codon_sequence):
			if codon in codon_to_position:
				codon_to_position[codon].append(i)
			else:
				codon_to_position['other'].append(i)

		# Get data
		variant_to_time = {}
		variant_to_v_max = {}
		variant_to_f_charged = {}
		variant_to_ribosome_position = {}
		variant_to_ribosome_per_mrna = {}
		variant_to_ribosome_termination = {}
		variant_to_argA_monomer = {}
		variant_to_argA_complex = {}
		variant_to_reaction = {}
		variant_to_arg = {}
		variant_to_amino_acid_saturation = {}
		variant_to_charging_rate = {}
		variant_to_ribosome_rate = {}
		variant_to_protein = {}

		variant_to_division_times = {}

		n_assume_feasible = 46
		cell_density = sim_data.constants.cell_density.asNumber(
			MASS_UNITS / VOLUME_UNITS)
		mmol_per_g_per_h = units.mmol / units.g / units.h

		for variant, seed in zip(variants, [8, 8]):
			cell_paths = self.ap.get_cells(variant=[variant], seed=[seed])

			# Time
			variant_to_time[variant] = read_stacked_columns(
				cell_paths, 'Main', 'time', remove_first=True).reshape(-1) / 60

			# Number of molecules
			(n_molecules,) = read_stacked_bulk_molecules(
				cell_paths, (molecules,), remove_first=True)

			# Max charging rate
			volume = 1e-15 * units.L * read_stacked_columns(
				cell_paths, 'Mass', 'cellVolume', remove_first=True).reshape(-1)
			c_synthetase = (1
				/ sim_data.constants.n_avogadro
				/ volume
				* n_molecules[:, molecules.index(synthetase)]
				).asNumber(units.umol / units.L)
			k_cat = (sim_data.relation.synthetase_to_k_cat[synthetase]
				).asNumber(1/units.s)
			variant_to_v_max[variant] = k_cat * c_synthetase

			# Charged fraction
			indexes = [molecules.index(trna) for trna in free_trnas]
			n_free_trnas = n_molecules[:, indexes].sum(axis=1)
			indexes = [molecules.index(trna) for trna in charged_trnas]
			n_charged_trnas = n_molecules[:, indexes].sum(axis=1)
			variant_to_f_charged[variant] = (n_charged_trnas
				/ (n_free_trnas + n_charged_trnas))

			# Ribosome position
			variant_to_division_times[variant] = []
			ribosomes = np.zeros((
				variant_to_time[variant].shape[0],
				len(codon_sequence),
				), dtype=np.int)

			for sim_dir in cell_paths:
				sim_out_dir = os.path.join(sim_dir, 'simOut')
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				# time_step = main_reader.readColumn('timeStepSec')
				time = main_reader.readColumn('time') / 60
				variant_to_division_times[variant].append(time[-1])

				trna_reader = TableReader(
					os.path.join(sim_out_dir, 'TrnaCharging'))
				ribosome_positions = trna_reader.readColumn(
					'ribosome_positions_argA')
				collision_removal = trna_reader.readColumn(
					'collision_removed_ribosomes_argA')
				ribosome_initiations = trna_reader.readColumn(
					'ribosome_initiation_argA')				

				for j in range(2, ribosome_positions.shape[0]):
					# Get index
					i = np.where(variant_to_time[variant] == time[j])[0][0]

					# Get current ribosome positions
					current = []
					for k, n in enumerate(ribosome_positions[j]):
						for l in range(n):
							current.append(k)
					original_current = copy.deepcopy(current)

					# No active ribosomes
					if len(current) == 0:
						continue

					# # Get current time step
					# time = time_step[j]

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
						k = previous.index(x)
						_ = previous.pop(k)

					# Get collision removals
					current_collision_removed = []
					for k, n in enumerate(collision_removal[j]):
						for l in range(n):
							current_collision_removed.append(k)

					current = np.array(current)
					previous = np.array(previous)

					# Newly initiated ribosomes
					n_initiations = np.logical_not(
						np.isnan(ribosome_initiations[j - 1])).sum()

					for k in range(n_initiations):
						_current = current[0]
						steps = _current + 1
						t = time / steps
						# ribosomes[j, 0:_current + 1] += t
						ribosomes[i, 0:_current + 1] += 1
						current = current[1:]

					for _previous in previous:

						# Skip terminating ribosomes
						if _previous == ribosome_positions.shape[1] - 1:
							continue

						elif np.any(current >= _previous):
							k = np.where(current >= _previous)[0][0]
							_current = current[k]
							steps = _current - _previous

							# Paused ribosomes
							if steps == 0:
								# ribosomes[j, _previous] += time
								ribosomes[i, _previous] += 1
								current = np.hstack((
									current[0:k], current[k + 1:]))

							# Processed ribosomes
							elif steps <= n_assume_feasible:
								t = time / steps
								# ribosomes[j, _previous + 1:_current + 1] += t
								ribosomes[i, _previous + 1:_current + 1] += 1
								current = np.hstack((
									current[0:k], current[k + 1:]))

							else:
								print('edge case 1: increase n_assume_feasible?')
								pass

						# Processing, then terminated ribosomes
						elif _previous + n_assume_feasible >= ribosome_positions.shape[1] - 1:
							_current = ribosomes.shape[1] - 1
							steps = _current - _previous
							t = time / steps
							# ribosomes[j, _previous + 1:] += t
							ribosomes[i, _previous + 1:] += 1

						else:
							print('edge case 2')
							pass

			variant_to_ribosome_position[variant] = ribosomes

			# Ribosomes per argA mRNA
			ribosome_positions = read_stacked_columns(
				cell_paths, 'TrnaCharging', 'ribosome_positions_argA',
				remove_first=True)[:, :-1]
			n_mrna = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_counts',
				remove_first=True)[:, mrna_index]
			mrna_present = n_mrna != 0

			ribosomes_per_mrna = np.zeros_like(n_mrna)
			ribosomes_per_mrna[mrna_present] = (
				ribosome_positions[mrna_present].sum(axis=1)
				/ n_mrna[mrna_present])
			variant_to_ribosome_per_mrna[variant] = ribosomes_per_mrna

			# Premature ribosome termination rate
			ribosome_terminations = read_stacked_columns(
				cell_paths, 'TrnaCharging', 'collision_removed_ribosomes_argA',
				remove_first=True)[:, :-1].sum(axis=1)
			time_step = units.s * read_stacked_columns(
				cell_paths, 'Main', 'timeStepSec', remove_first=True,
				).reshape(-1)
			ribosome_termination_rate = (1
				/ time_step
				* ribosome_terminations
				).asNumber(1/units.s)
			variant_to_ribosome_termination[variant] = ribosome_termination_rate

			# argA enzyme monomer and complex counts
			variant_to_argA_monomer[variant] = n_molecules\
				[:, molecules.index(enzyme_monomer)]
			variant_to_argA_complex[variant] = n_molecules\
				[:, molecules.index(enzyme_complex)]

			# Flux through argA reaction
			cell_mass = read_stacked_columns(
				cell_paths, 'Mass', 'cellMass', remove_first=True).reshape(-1)
			dry_mass = read_stacked_columns(
				cell_paths, 'Mass', 'dryMass', remove_first=True).reshape(-1)
			raw_reaction_fluxes = read_stacked_columns(
				cell_paths, 'FBAResults', 'reactionFluxes', remove_first=True)
			coefficient = dry_mass / cell_mass * cell_density

			sim_out_dir = os.path.join(cell_paths[0], 'simOut')
			fba_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
			reaction_ids = fba_reader.readAttribute('reactionIDs')
			reaction_index = reaction_ids.index(reaction)
			reaction_flux = (
				COUNTS_UNITS
				/ MASS_UNITS
				/ TIME_UNITS
				* raw_reaction_fluxes[:, reaction_index]
				/ coefficient).asNumber(mmol_per_g_per_h)
			variant_to_reaction[variant] = reaction_flux

			# Arginine concentration
			c_arg = (1
				/ sim_data.constants.n_avogadro
				/ volume
				* n_molecules[:, molecules.index(amino_acid)]
				).asNumber(units.umol / units.L)
			variant_to_arg[variant] = c_arg

			# Arginine saturation fraction
			K_A = (sim_data.relation.synthetase_to_K_A[synthetase]
				).asNumber(units.umol/units.L)
			variant_to_amino_acid_saturation[variant] = (c_arg
				/ (K_A + c_arg))

			# Charging rate
			charging_events = read_stacked_columns(
				cell_paths, 'TrnaCharging', 'charging_events',
				remove_first=True)[:, trna_indexes].sum(axis=1)
			charging_rate = (1
				/ sim_data.constants.n_avogadro
				/ volume
				/ time_step
				* charging_events
				).asNumber(units.umol / units.L / units.s)
			variant_to_charging_rate[variant] = charging_rate

			# Ribosome elongation rate
			variant_to_ribosome_rate[variant] = read_stacked_columns(
				cell_paths, 'RibosomeData', 'effectiveElongationRate',
				remove_first=True).reshape(-1)

			# Protein mass
			variant_to_protein[variant] = read_stacked_columns(
				cell_paths, 'Mass', 'proteinMass',
				remove_first=True)

		################################################################
		# Plot
		fig, axes = plt.subplots(
			13, 2, figsize=(7.5, 7), sharey='row', sharex='col')

		color_control = '#bdbdbd'
		color_experiment = '#6b6ecf'

		color_non_coi = '#e7e7e7'# '#d9d9d9'
		colors = [
			'#6b6ecf', # blue
			'#b5cf6b', # green
			'#e7ba52', # yellow
			'#ce6dbd', # pink
			]
		codons = ['CGU', 'CGG', 'CGC', 'AGG']
		codon_to_colors = dict(zip(codons, colors))

		color_highlight_red = '#d6616b'

		fontsize_legend = 7
		fontsize_label = 8
		fontsize_heading = 11
		linewidth = 1

		window_time = 60 * 10 # seconds
		window_size = window_time // 2 # 2 second time steps
		window = np.ones(window_size)

		def plot(axes, variant_to_data, highlight_zero=False, y=0, smoothen=False):
			# Control
			ax = axes[0]
			ax.plot(
				variant_to_time[variant_control],
				variant_to_data[variant_control],
				linewidth=linewidth,
				color=color_control,
				)

			# Experiment
			ax = axes[1]
			if smoothen:
				smoothened = np.convolve(
					variant_to_data[variant_experiment], window, mode='same') / window_size
				ax.plot(
					variant_to_time[variant_experiment][window_size:-window_size],
					smoothened[window_size:-window_size],
					linewidth=linewidth,
					color=color_experiment,
					)
			else:
				ax.plot(
					variant_to_time[variant_experiment],
					variant_to_data[variant_experiment],
					linewidth=linewidth,
					color=color_experiment,
					)

			if highlight_zero:
				where_zero = np.where(variant_to_data[variant_experiment] == 0)[0]
				for i in range(where_zero.shape[0]):
					continuous_array = np.arange(where_zero[i], where_zero[-1] + 1)
					if np.all(where_zero[i:] == continuous_array):
						ax.plot(
							[variant_to_time[variant_experiment][where_zero[0]],
							variant_to_time[variant_experiment][where_zero[-1]]],
							[y, y],
							color_highlight_red,
							)
						break
			return

		label_prop = {
			'fontsize': fontsize_label,
			'rotation': 'horizontal',
			'labelpad': 60.,
		}

		scatter_prop = {
			'marker': '.',
			's': 1,
			'linewidths': 0,
		}

		# Max charging rate
		row = 0
		axes[row, 0].set_ylabel(
			'Maximal Aminoacylation\nRate (uM/s)', **label_prop)
		plot(axes[row, :], variant_to_v_max)
		print('Max charging rate')
		print(f'\tOptimized:\t{variant_to_v_max[0].mean():.1f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		print(f'\tMeasured:\t{variant_to_v_max[6][:i].mean():.1f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][12])[0][0]
		j = np.where(variant_to_time[6] == variant_to_division_times[6][15])[0][0]
		print(f'\t\tGen 12 to 15: {variant_to_v_max[6][i:j].mean():.1f}')

		# Text labels
		y = 250
		for ax, variant in zip(axes[row, :], variants):
			for i, division_time in enumerate(variant_to_division_times[variant]):
				if i == 0:
					tau = division_time
				else:
					tau = division_time - variant_to_division_times[variant][i -1]
				x = division_time - (tau / 2)

				if i % 2 == 0:
					if i == 0:
						text = f'Gen\n0'
					else:
						text = f'{i}'
					ax.text(x, y, text, fontsize=fontsize_legend,
						va='bottom', ha='center')

		# Charged fraction
		row += 1
		axes[row, 0].set_ylabel(
			'Aminoacylated Fraction', **label_prop)
		plot(axes[row, :], variant_to_f_charged)
		print('Aminoacylated fraction')
		print(f'\tOptimized:\t{variant_to_f_charged[0].mean():.3f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][13])[0][0]
		print(f'\tMeasured:\t{variant_to_f_charged[6][:i].mean():.3f}')

		################################################################
		# Ribosome position
		row += 1
		axes[row, 0].set_ylabel(
			'Ribosome A Site Position\n(sequence)', **label_prop)

		for ax, variant in zip(axes[row, :], variants):
			time = variant_to_time[variant]
			x, y = np.where(variant_to_ribosome_position[variant])

			# Non arginine codons
			mask = [_y in codon_to_position['other'] for _y in y]
			ax.scatter(
				time[x[mask]], y[mask],
				color=color_non_coi,
				**scatter_prop,
				rasterized=True,
				)

			# Arginine codons
			for codon in codons:
				color = codon_to_colors[codon]
				mask = [_y in codon_to_position[codon] for _y in y]
				ax.scatter(
					time[x[mask]], y[mask],
					color=color,
					**scatter_prop,
					rasterized=True,
					)

			# Format
			ax.set_yticks([0, len(codon_sequence) - 1])

			# Indicate CGG positions
			if variant == variant_experiment:
				to_highlight = []
				for position in codon_to_position['CGG']:
					if 100 < position < 200:
						to_highlight.append(position)
				for position in to_highlight:
					ax.plot([1800, 2000], [position, position], color='k', linewidth=0.5)

		################################################################
		# Ribosomes on argA mrnas
		row += 1
		axes[row, 0].set_ylabel(
			'Ribosome Accumulation\n(ribosomes/mRNA)', **label_prop)
		plot(axes[row, :], variant_to_ribosome_per_mrna)
		for ax in axes[row, :]:
			ax.set_yticks([0, 40])
		print('Ribosome accumulation')
		print(f'\tOptimized:\t{variant_to_ribosome_per_mrna[0].mean():.1f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		print(f'\tMeasured:\t{variant_to_ribosome_per_mrna[6][:i].mean():.1f}')
		print(f'\t\tMax: {np.max(variant_to_ribosome_per_mrna[6][:i])}')

		# Premature ribosome termination rate
		row += 1
		axes[row, 0].set_ylabel(
			'Premature Ribosome\nTermination Rate (events/s)', **label_prop)
		plot(axes[row, :], variant_to_ribosome_termination)
		for ax in axes[row, :]:
			ax.set_yticks([0, 18])
		print('Premature ribosome termination')
		print(f'\tOptimized:\t{variant_to_ribosome_termination[0].mean():.3f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		print(f'\tMeasured:\t{variant_to_ribosome_termination[6][:i].mean():.3f}')
		print(f'\t\tMax: {np.max(variant_to_ribosome_termination[6][:i]):.3f}')


		################################################################
		# argA monomer counts
		row += 1
		axes[row, 0].set_ylabel('Monomer\n(molecules/cell)', **label_prop)
		plot(axes[row, :], variant_to_argA_monomer)
		for ax in axes[row, :]:
			ax.set_yticks([0, 6, max(variant_to_argA_monomer[0])])
		print('ArgA monomers')
		print(f'\tOptimized:\t{variant_to_argA_monomer[0].mean():.2f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][7])[0][0]
		print(f'\tMeasured:\t{variant_to_argA_monomer[6][:i].mean():.2f}')

		# argA complex counts
		row += 1
		axes[row, 0].set_ylabel('Hexamer\n(molecules/cell)', **label_prop)
		plot(axes[row, :], variant_to_argA_complex, highlight_zero=True, y=200)
		print('ArgA hexamers')
		print(f'\tOptimized:\t{variant_to_argA_complex[0].mean():.2f}')
		i_divisions = [np.where(variant_to_time[0] == x)[0][0] for x in variant_to_division_times[0]]
		fold_changes = []
		for i, i_division in enumerate(i_divisions):
			if i == 0:
				i_start = 0
			else:
				i_start = i_divisions[i - 1] + 1
			fold_changes.append(variant_to_argA_complex[0][i_division] / variant_to_argA_complex[0][i_start])
		print(f'\t\tAvg cell cycle increase: {np.mean(fold_changes):.1f} fold')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		print(f'\tMeasured:\t{variant_to_argA_complex[6][:i].mean():.2f}')

		# argA reaction
		row += 1
		axes[row, 0].set_ylabel(
			'N-Acetyl Transfer\nReaction (mmol/gDCW/h)', **label_prop)
		plot(axes[row, :], variant_to_reaction,
			highlight_zero=True, y=1e2, smoothen=True)
		for ax in axes[row]:
			ax.set_yscale('log')
			ax.set_yticks([1e-2, 1e2])
			ax.set_ylim([1e-2, 2e2])
		print('N-acetyl transfer reaction')
		print(f'\tOptimized:\t{variant_to_reaction[0].mean():.3f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		print(f'\tMeasured:\t{variant_to_reaction[6][:i].mean():.3f}')

		# Arginine
		row += 1
		axes[row, 0].set_ylabel('Arginine (uM)', **label_prop)
		plot(axes[row, :], variant_to_arg)
		y_avg = round(variant_to_arg[0].mean() / 10) * 10
		for ax in axes[row, :]:
			ax.set_yticks([0, y_avg, 1000])
		print('Arginine')
		print(f'\tOptimized:\t{variant_to_arg[0].mean():.1f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		j = np.where(variant_to_time[6] == variant_to_division_times[6][13])[0][0]
		print(f'\tMeasured:\t{variant_to_arg[6][i:j].mean():.1f}')

		# Arginine saturation fraction
		row += 1
		axes[row, 0].set_ylabel('ArgRS Fractional Saturation\nfor Arginine', **label_prop)
		plot(axes[row, :], variant_to_amino_acid_saturation)
		for ax in axes[row, :]:
			ax.set_yticks([0, 0.5, 1])
			ax.set_ylim([0, 1])
		print('ArgRS fractional saturation for arginine')
		print(f'\tOptimized:\t{variant_to_amino_acid_saturation[0].mean():.3f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		j = np.where(variant_to_time[6] == variant_to_division_times[6][13])[0][0]
		print(f'\tMeasured:\t{variant_to_amino_acid_saturation[6][i:j].mean():.3f}')

		# Charging rate
		row += 1
		axes[row, 0].set_ylabel('Aminoacylation Rate\n(uM/s)', **label_prop)
		plot(axes[row, :], variant_to_charging_rate, highlight_zero=True, y=20)
		print('Aminoacylation rate')
		print(f'\tOptimized:\t{variant_to_charging_rate[0].mean():.1f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		j = np.where(variant_to_time[6] == variant_to_division_times[6][13])[0][0]
		print(f'\tMeasured:\t{variant_to_charging_rate[6][:i].mean():.2f} (Gen 0 through 8)')
		print(f'\tMeasured:\t{variant_to_charging_rate[6][i:j].mean():.2f} (Gen 9 through 13)')

		# Ribosome rate
		row += 1
		axes[row, 0].set_ylabel(
			'Ribosome Elongation Rate\n(amino acids/s/ribosome)', **label_prop)
		plot(axes[row, :], variant_to_ribosome_rate, highlight_zero=True, y=20)
		for ax in axes[row, :]:
			ax.set_yticks([0, 10, 20])
		print('Ribosome elongation rate')
		print(f'\tOptimized:\t{variant_to_ribosome_rate[0].mean():.1f}')
		i = np.where(variant_to_time[6] == variant_to_division_times[6][8])[0][0]
		j = np.where(variant_to_time[6] == variant_to_division_times[6][13])[0][0]
		print(f'\tMeasured:\t{variant_to_ribosome_rate[6][:i].mean():.2f} (Gen 0 through 8)')
		print(f'\tMeasured:\t{variant_to_ribosome_rate[6][i:j].mean():.2f} (Gen 9 through 13)')

		# Protein mass
		row += 1
		axes[row, 0].set_ylabel('Protein Mass (fg)', **label_prop)
		plot(axes[row, :], variant_to_protein)
		mass_initial = [variant_to_protein[0][0][0]]
		mass_final = []
		for x in variant_to_division_times[0]:
			i = np.where(variant_to_time[0] == x)[0][0]

			mass_final.append(variant_to_protein[0][i][0])
			if i + 1 < len(variant_to_protein[0]):
				mass_initial.append(variant_to_protein[0][i + 1][0])
		y_min = round(np.mean(mass_initial) / 10) * 10
		y_max = round(np.mean(mass_final) / 10) * 10
		for ax in axes[row, :]:
			ax.set_yticks([0, y_min, y_max])

		################################################################

		# Division times
		for ax in axes[:, 0]:
			for x in variant_to_division_times[0]:
				ax.axvline(x, color='k', linewidth=0.2, linestyle='--')
		for ax in axes[:, 1]:
			for x in variant_to_division_times[6]:
				ax.axvline(x, color='k', linewidth=0.2, linestyle='--')

		# Format
		for ax in axes[-1, :]:
			ax.set_xlabel('Time (min)', fontsize=fontsize_label)

		y_portion = 0.08
		for ax in axes.reshape(-1):
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(
				axis='both', which='major', labelsize=fontsize_legend)
			ax.margins(0.01)

			y_max = ax.get_ylim()[1]
			y_min = -(y_max * y_portion / (1 - y_portion))
			ax.set_ylim([y_min, y_max])

			# Premature ribosome termination
			if ax in axes[4, :]:
				ax.set_ylim([ax.get_ylim()[0], 18])

		plt.tight_layout(h_pad=0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
