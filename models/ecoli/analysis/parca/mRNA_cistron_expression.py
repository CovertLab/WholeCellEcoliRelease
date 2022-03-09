"""
Plots the differences between the relative expression levels of each mRNA
cistron that is expected from RNAseq data and the actual expression levels of
each mRNA cistron after accounting for operon structure.
"""

import pickle

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure


CONDITION = 'basal'
NUMERICAL_ZERO = 1e-20
LABEL_BOUNDARY = 10
PLOT_OUTLIER_TUS = False
SUPPLEMENTAL_PLOT_EXT = '_tu_structures'


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load data from sim_data
		all_cistron_ids = sim_data.process.transcription.cistron_data['id']
		all_rna_ids = sim_data.process.transcription.rna_data['id']
		cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		cistron_is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		is_ribosomal_protein = sim_data.process.transcription.cistron_data[
			'is_ribosomal_protein']
		is_rnap = sim_data.process.transcription.cistron_data['is_RNAP']

		# Get mask for genes that are
		# i) mRNAs
		# ii) Does not encode for ribosomal proteins or RNAPs
		# iii) not manually overexpressed
		# to isolate the effects of operons on expression levels.
		is_adjusted = np.zeros_like(cistron_is_mRNA, dtype=bool)
		for adjusted_cistron_id in sim_data.adjustments.rna_expression_adjustments.keys():
			# Include cistrons whose expression is adjusted because they belong
			# to the same TU as the cistron that is bumped up
			adjusted_rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(
				adjusted_cistron_id)
			adjusted_cistron_indexes = []
			for adjusted_rna_index in adjusted_rna_indexes:
				adjusted_cistron_indexes.extend(
					sim_data.process.transcription.rna_id_to_cistron_indexes(
						all_rna_ids[adjusted_rna_index]))

			is_adjusted[adjusted_cistron_indexes] = True

		mask = cistron_is_mRNA & ~is_ribosomal_protein & ~is_rnap & ~is_adjusted

		# Get boolean array of relevant cistrons that belong to at least one
		# polycistronic transcript
		polycistronic_cistron_indexes = []
		for rna_id in sim_data.process.transcription.rna_data['id']:
			cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes(rna_id)
			if len(cistron_indexes) > 1:
				polycistronic_cistron_indexes.extend(cistron_indexes)
		is_polycistronic = np.zeros(len(all_cistron_ids), bool)
		if len(polycistronic_cistron_indexes) > 0:
			is_polycistronic[np.array(
				list(set(polycistronic_cistron_indexes)))] = True
		is_polycistronic = is_polycistronic[mask]

		# Get expression
		expected_mRNA_cistron_exp = sim_data.process.transcription.cistron_expression[CONDITION]
		rna_exp = sim_data.process.transcription.rna_expression[CONDITION]
		actual_mRNA_cistron_exp = cistron_tu_mapping_matrix.dot(rna_exp)

		# Normalize expression
		expected_mRNA_cistron_exp /= expected_mRNA_cistron_exp[mask].sum()
		actual_mRNA_cistron_exp /= actual_mRNA_cistron_exp[mask].sum()

		# Find cistrons with more than a 10-fold difference between actual vs
		# expected expression
		exp_diff = np.log10(actual_mRNA_cistron_exp + NUMERICAL_ZERO) - np.log10(expected_mRNA_cistron_exp + NUMERICAL_ZERO)
		largest_diff_indexes = np.where(np.abs(exp_diff) > np.log10(LABEL_BOUNDARY))[0]
		largest_diff_indexes = [index for index in largest_diff_indexes if mask[index]]

		plt.figure(figsize=(12, 12))
		ls = np.logspace(-8, -1, 8)
		plt.plot(ls, LABEL_BOUNDARY * ls, c='#dddddd', ls='--')
		plt.plot(ls, 1 / LABEL_BOUNDARY * ls, c='#dddddd', ls='--')

		plt.scatter(
			expected_mRNA_cistron_exp[mask], actual_mRNA_cistron_exp[mask],
			c='#dddddd', s=2, label='monocistronic')

		# Highlight genes that are polycistronic
		plt.scatter(
			expected_mRNA_cistron_exp[mask][is_polycistronic],
			actual_mRNA_cistron_exp[mask][is_polycistronic],
			c='#333333', s=4, label='polycistronic')

		# Label cistrons that have more than a 10-fold difference
		for index in largest_diff_indexes:
			plt.text(
				expected_mRNA_cistron_exp[index],
				actual_mRNA_cistron_exp[index] * 1.05,
				all_cistron_ids[index],
				ha='center', va='bottom', fontsize=5)

		plt.xlabel('Original expression expected from RNAseq')
		plt.ylabel('Actual expression after applying operon structure')
		plt.xlim([1e-8, 1e-1])
		plt.ylim([1e-8, 1e-1])
		plt.xscale('log')
		plt.yscale('log')
		plt.legend()

		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')

		if PLOT_OUTLIER_TUS:
			# Load data from sim_data
			cistron_id_to_gene_id = {
				cistron['id']: cistron['gene_id']
				for cistron in sim_data.process.transcription.cistron_data
				}
			tu_ids = [x[:-3] for x in sim_data.process.transcription.rna_data['id']]

			# Divide mapping matrix into individual operons
			cistron_indexes, tu_indexes = cistron_tu_mapping_matrix.nonzero()

			# Sort operons by size (number of cistrons)
			all_operons = sorted(
				sim_data.process.transcription.operons,
				key=lambda operon: len(operon[0]), reverse=True)

			# Create map of cistron index to operon index
			cistron_to_operon = {}
			for (operon_index, (operon_cistron_indexes, _)
				 ) in enumerate(all_operons):
				for cistron_index in operon_cistron_indexes:
					cistron_to_operon[cistron_index] = operon_index


			# Create dictionary containing expected and actual expression data for all cistrons
			# in operons with at least one cistron with more than a 10-fold difference between actual vs
			# expected expression
			diff_expression_operon_dict = {}
			bad_cistrons = []
			for i, cistron_diff_index in enumerate(largest_diff_indexes):
				try:
					diff_operon_index = cistron_to_operon[cistron_diff_index]
					operon_cistron_indexes, operon_tu_indexes = all_operons[diff_operon_index]
					cistron_labels = []
					expected_exp = []
					actual_exp = []
					for cistron_index in operon_cistron_indexes:
						cistron_labels.append(
							cistron_id_to_gene_id[all_cistron_ids[cistron_index]])
						expected_exp.append(
							expected_mRNA_cistron_exp[cistron_index])
						actual_exp.append(actual_mRNA_cistron_exp[cistron_index])
					diff_expression_operon_dict[i] = (
						operon_cistron_indexes,
						operon_tu_indexes,
						cistron_labels,
						expected_exp,
						actual_exp)
				except Exception as e:
					print(f'Scanning cistrons: {e}')
					bad_cistrons.append(cistron_diff_index)

			fig = plt.figure(figsize=(8, 8))
			fig.set_size_inches(50, 50)

			gs = gridspec.GridSpec(
				len(diff_expression_operon_dict)*2 // 10 + 1, 10, figure=fig)
			for i in range(len(diff_expression_operon_dict)):
				operon_cistron_indexes, operon_tu_indexes, cistron_labels, expected_exp, actual_exp = diff_expression_operon_dict[i]
				# Plot actual and expected expression data for all cistrons in
				# operon
				ax1 = fig.add_subplot(gs[i * 2 // 10, i * 2 % 10])
				x = np.arange(len(cistron_labels))
				width = 0.35
				ax1.bar(x - width / 2, expected_exp, width, log=True, label='expected')
				ax1.bar(x + width / 2, actual_exp, width, log=True, label='actual')
				ax1.set_xticks(x)
				ax1.set_xticklabels(cistron_labels)
				ax1.legend()
				# Get TU structure array of this operon
				tu_structure = np.full(
					(len(operon_cistron_indexes), len(operon_tu_indexes)), 255)
				mask = np.isin(cistron_indexes, operon_cistron_indexes)
				for (x, y, z) in zip(
					cistron_indexes[mask],
					tu_indexes[mask],
						cistron_tu_mapping_matrix.data[mask]):
					tu_structure[
						np.where(operon_cistron_indexes == x)[0][0],
						np.where(operon_tu_indexes == y)[0][0]] = 0

				# Display TU structure as binary heatmap
				ax2 = fig.add_subplot(gs[(i * 2 + 1) // 10, (i * 2 + 1) % 10])
				ax2.imshow(tu_structure, cmap='gray')
				ax2.set_xticks(np.arange(len(operon_tu_indexes)))
				ax2.set_yticks(np.arange(len(operon_cistron_indexes)))
				ax2.set_xticklabels(
					[tu_ids[i] for i in operon_tu_indexes],
					rotation=90, fontsize=6)
				ax2.set_yticklabels(
					[cistron_id_to_gene_id[all_cistron_ids[i]] for i in operon_cistron_indexes],
					fontsize=6)
				ax2.spines['top'].set_visible(False)
				ax2.spines['right'].set_visible(False)
				ax2.spines['bottom'].set_visible(False)
				ax2.spines['left'].set_visible(False)
				ax2.tick_params(axis='both', bottom=False, left=False, pad=1)

			plt.tight_layout()
			exportFigure(plt, plot_out_dir, plot_out_filename + SUPPLEMENTAL_PLOT_EXT, metadata)
			plt.close('all')


if __name__ == "__main__":
	Plot().cli()
