"""
Plots trna charging specs.
"""

import pickle
import os

from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


trna_to_label = {
	'RNA0-300[c]': 'valZ',
	'RNA0-301[c]': 'lysY',
	'RNA0-302[c]': 'lysZ',
	'RNA0-303[c]': 'lysQ',
	'RNA0-304[c]': 'asnW',
	'RNA0-305[c]': 'ileY',
	'RNA0-306[c]': 'metV',
	}

amino_acid_to_label = {
	'L-ALPHA-ALANINE[c]': 'Alanine',
	'ARG[c]': 'Arginine',
	'ASN[c]': 'Asparagine',
	'L-ASPARTATE[c]': 'Aspartate',
	'CYS[c]': 'Cysteine',
	'GLT[c]': 'Glutamate',
	'GLN[c]': 'Glutamine',
	'GLY[c]': 'Glycine',
	'HIS[c]': 'Histidine',
	'ILE[c]': 'Isoleucine',
	'LEU[c]': 'Leucine',
	'LYS[c]': 'Lysine',
	'MET[c]': 'Methionine',
	'PHE[c]': 'Phenylalanine',
	'PRO[c]': 'Proline',
	'SER[c]': 'Serine',
	'THR[c]': 'Threonine',
	'TRP[c]': 'Tryptophan',
	'TYR[c]': 'Tyrosine',
	'VAL[c]': 'Valine',
	}

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		with open(self.ap.get_variant_kb(variants[0]), 'rb') as f:
			variant_sim_data = pickle.load(f)

		if not hasattr(variant_sim_data, 'trna_synthetase_kinetics_variant'):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics')
			return

		# Analyze variants 0 (control) and 1 (previous)
		if not (0 in variants and 1 in variants):
			print('This analysis script is designed for the variant:'
				' trna_synthetase_kinetics variant indexes 0 and 1.')
			return

		# Describe molecules
		amino_acid = 'L-ALPHA-ALANINE[c]'
		free_trna = 'alaT-tRNA[c]'

		# trna indexes
		rna_data = sim_data.process.transcription.rna_data
		all_trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		trna_index = all_trnas.index(free_trna)
		charged_trna = sim_data.process.transcription.charged_trna_names[trna_index]

		# Alanine codons
		ala_codons = sim_data.relation.amino_acid_to_codons[amino_acid]
		ala_codon_indexes = [sim_data.relation.codons.index(codon)
			for codon in ala_codons]

		# alaT-trna codons
		alaT_codons = []
		for codon in ala_codons:
			if charged_trna in sim_data.relation.codon_to_trnas[codon]:
				alaT_codons.append(codon)
		alaT_codon_indexes = [ala_codons.index(codon) for codon in alaT_codons]

		# alaT-trna-codon interactions
		trna_codon_indexes = [sim_data.relation.trna_codon_pairs
			.index(f'{free_trna}_{codon}') for codon in alaT_codons]

		# Plot
		def make_plot(ax, y, color, shade_color, label='', skip_first=False):
			i = 0
			if skip_first:
				i = 1

			y_mean = y[:, i:].mean(axis=0)
			y_std = y[:, i:].std(axis=0)
			ax.plot(relative_time[i:], y_mean, color=color, label=label)
			ax.fill_between(relative_time[i:], y_mean + y_std, y_mean - y_std,
				color=shade_color, alpha=0.4)
			return

		def make_blank(ax):
			ax.text(0.5, 0.5, 'Not\nRepresented', ha='center', va='center',
				transform=ax.transAxes, fontsize=legend_fontsize)
			ax.set_xticks([])
			ax.set_yticks([])
			return

		fig, axes = plt.subplots(2, 6, figsize=(8, 2.75))
		label_fontsize = 7
		legend_fontsize = 7

		colors = ['#5254a3', '#8ca252', '#bd9e39', '#a55194']
		shading_colors = ['#9c9ede', '#cedb9c', '#e7cb94', '#de9ed6']
		color_control = '#969696'
		color_control_shading = '#d9d9d9'
		free_trna_label = f'{trna_to_label.get(free_trna, free_trna[:4])}-tRNA'
		amino_acid_label = amino_acid_to_label[amino_acid]

		################################################################

		variants = [0, 1]
		for variant in variants:
			# Load modified variant sim_data
			## Consider calculating variant difference and only loading
			## sim_data once above for better performance.
			with open(self.ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)

			cell_paths = self.ap.get_cells(variant=[variant])

			relative_time = np.linspace(0, 1, 1000)
			free_trna_array = []
			charged_trna_array = []
			total_trna_array = []
			charged_fraction_array = []
			charging_rates_array = []
			trna_codon_interaction_array = []
			alaT_codons_read_array = []
			ala_codons_read_array = []
			aas_read_array = []
			cleaved_array = []

			for sim_dir in cell_paths:
				sim_out_dir = os.path.join(sim_dir, 'simOut')

				# Listeners used
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				charging_reader = TableReader(os.path.join(sim_out_dir, 'TrnaCharging'))

				# Load data
				time = main_reader.readColumn('time')
				time -= time[0]
				time /= time[-1]
				time_step = main_reader.readColumn('timeStepSec')

				(n_free_trna, n_charged_trna) = read_bulk_molecule_counts(
					sim_out_dir, ([free_trna], [charged_trna]))

				charging_events = charging_reader.readColumn(
					'charging_events')[:, trna_index]
				charging_rate = charging_events / time_step

				trna_codon_interaction = charging_reader.readColumn(
					'codons_to_trnas_counter')[:, trna_codon_indexes] / time_step[:, None]

				all_codons_read = charging_reader.readColumn('codons_read')
				ala_codons_read = all_codons_read[:, ala_codon_indexes] / time_step[:, None]
				alaT_codons_read = ala_codons_read[:, alaT_codon_indexes] / time_step[:, None]

				cleaved = charging_reader.readColumn('cleaved') / time_step

				# Normalize
				f = interp1d(time, n_free_trna, copy=False)
				free_trna_array.append(f(relative_time))

				f = interp1d(time, n_charged_trna, copy=False)
				charged_trna_array.append(f(relative_time))

				total_trna = n_free_trna + n_charged_trna
				f = interp1d(time, total_trna, copy=False)
				total_trna_array.append(f(relative_time))

				charged_fraction = n_charged_trna / total_trna
				f = interp1d(time, charged_fraction, copy=False)
				charged_fraction_array.append(f(relative_time))

				f = interp1d(time, charging_rate, copy=False)
				charging_rates_array.append(f(relative_time))

				f = interp1d(time, trna_codon_interaction, axis=0, copy=False)
				trna_codon_interaction_array.append(f(relative_time))

				f = interp1d(time, alaT_codons_read, axis=0, copy=False)
				alaT_codons_read_array.append(f(relative_time))

				f = interp1d(time, ala_codons_read, axis=0, copy=False)
				ala_codons_read_array.append(f(relative_time))

				aas_read = ala_codons_read.sum(axis=1)
				f = interp1d(time, aas_read, axis=0, copy=False)
				aas_read_array.append(f(relative_time))

				f = interp1d(time, cleaved, axis=0, copy=False)
				cleaved_array.append(f(relative_time))

			free_trna_array = np.array(free_trna_array)
			charged_trna_array = np.array(charged_trna_array)
			total_trna_array = np.array(total_trna_array)
			charged_fraction_array = np.array(charged_fraction_array)
			charging_rates_array = np.array(charging_rates_array)
			trna_codon_interaction_array = np.array(trna_codon_interaction_array)
			alaT_codons_read_array = np.array(alaT_codons_read_array)
			ala_codons_read_array = np.array(ala_codons_read_array)
			aas_read_array = np.array(aas_read_array)
			cleaved_array = np.array(cleaved_array)


			# Plot
			if variant == 1:

				# Previous
				ax = axes[0, 0]
				ax.set_title(f'{free_trna_label}\nAminoacylation', fontsize=label_fontsize)
				make_blank(ax)

				ax = axes[0, 1]
				ax.set_title(f'{free_trna_label}\nAbundance', fontsize=label_fontsize)
				ax.set_ylabel('Number of Molecules', fontsize=label_fontsize)
				make_plot(ax, total_trna_array, color_control, color_control_shading)
				ax.set_yticks([3e3, 7e3])

				ax = axes[0, 2]
				ax.set_title(f'{free_trna_label}\nAminoacylated Fraction', fontsize=label_fontsize)
				make_blank(ax)

				ax = axes[0, 3]
				ax.set_title(f'{free_trna_label}-Codon\nInteraction', fontsize=label_fontsize)
				make_blank(ax)

				ax = axes[0, 4]
				ax.set_title('Polypeptide\nElongation', fontsize=label_fontsize)
				ax.set_ylabel(f'{amino_acid_label}s / s', fontsize=label_fontsize)
				make_plot(ax, aas_read_array, color_control, color_control_shading, skip_first=True)
				ax.set_yticks([2e4, 5e4])

				ax = axes[0, 5]
				ax.set_title(f'N-terminal\nMethionine Cleavage', fontsize=label_fontsize)
				make_blank(ax)

			elif variant == 0:

				################################################################
				# Current
				ax = axes[1, 0]
				ax.set_ylabel('Events / s', fontsize=label_fontsize)
				make_plot(ax, charging_rates_array, colors[0], shading_colors[0], skip_first=True)
				ax.set_yticks([3e3, 9e3])
				# ax.set_ylim([3e3, 9e3])

				ax = axes[1, 1]
				ax.set_ylabel('Number of Molecules', fontsize=label_fontsize)
				make_plot(ax, charged_trna_array, colors[0], shading_colors[0])
				make_plot(ax, free_trna_array, color_control, color_control_shading)
				ax.set_yticks([0e3, 6e3])
				ax.text(1, 6e3, 'Aminoacylated', color=colors[0], ha='right', va='center', fontsize=legend_fontsize)
				ax.text(1, 1e3, 'Unaminoacylated', color=color_control, ha='right', va='center', fontsize=legend_fontsize)

				ax = axes[1, 2]
				ax.set_ylabel('Fraction', fontsize=label_fontsize)
				make_plot(ax, charged_fraction_array, colors[0], shading_colors[0])
				ax.set_yticks([0, 1])

				ax = axes[1, 3]
				ax.set_ylabel('Events / s', fontsize=label_fontsize)
				for i, codon in enumerate(alaT_codons):
					color = colors[ala_codons.index(codon)]
					shading_color = shading_colors[ala_codons.index(codon)]
					make_plot(ax, trna_codon_interaction_array[:, :, i], color,
						shading_color, skip_first=True)
				ax.set_yticks([2e3, 5e3])
				ax.text(1, 3e3, 'GCA', color=colors[0], ha='right', va='center', fontsize=legend_fontsize)
				ax.text(1, 4e3, 'GCG', color=colors[2], ha='right', va='center', fontsize=legend_fontsize)

				ax = axes[1, 4]
				ax.set_ylabel('Codons / s', fontsize=label_fontsize)
				for i, codon in enumerate(ala_codons):
					make_plot(ax, ala_codons_read_array[:, :, i], colors[i],
						shading_colors[i], skip_first=True)

					if codon == 'GCG':
						va = 'bottom'
					elif codon == 'GCA':
						va = 'top'
					else:
						va = 'center'
					ax.text(1.02, ala_codons_read_array[:, :, i][:, -1].mean(), codon,
						color=colors[i], ha='left', va=va, fontsize=legend_fontsize)
				ax.set_yticks([4e3, 1.8e4])

				ax = axes[1, 5]
				ax.set_ylabel('Events / s', fontsize=label_fontsize)
				make_plot(ax, cleaved_array, colors[0], shading_colors[0], skip_first=True)
				ax.set_yticks([4e2, 1e3])

		################################################################
		# Format
		for ax in axes.reshape(-1):
			ax.set_xticks([0, 1])

		for ax in [axes[1, 0], axes[0, 1], axes[1, 1], axes[1, 3], axes[0, 4], axes[1, 4], axes[1, 5]]:
			ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
			ax.yaxis.offsetText.set_fontsize(legend_fontsize)

		axes[1, 2].set_xlabel('Cell Cycle Progress', fontsize=label_fontsize)

		for ax in axes.reshape(-1):
			ax.tick_params(axis='both', which='major', labelsize=legend_fontsize)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		plt.tight_layout(w_pad=0.5)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
