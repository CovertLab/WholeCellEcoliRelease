"""
Plots the differences between the relative expression levels of each mRNA
cistron that is expected from RNAseq data and the actual expression levels of
each mRNA cistron after fitting RNA expression. These two values may be
different for one or more of these three reasons:

1) Adjustment of expression levels for genes encoding for RNAP/ribosome subunits
and the following normalization of expression levels
2) Manual adjustment of expression levels for specific genes listed in
reconstruction/flat/adjustments/rna_expression_adjustments.tsv.
3) Transcription unit structure forcing stoichiometric transcription of certain
cistrons.
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure


CONDITION = 'basal'
NUMERICAL_ZERO = 1e-20
LABEL_BOUNDARY = 10


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load data from sim_data
		cistron_ids = sim_data.process.transcription.cistron_data['id']
		cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		cistron_is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		mRNA_cistron_ids = cistron_ids[cistron_is_mRNA]
		mRNA_cistron_is_ribosomal_protein = sim_data.process.transcription.cistron_data['is_ribosomal_protein'][cistron_is_mRNA]
		mRNA_cistron_is_rnap = sim_data.process.transcription.cistron_data['is_RNAP'][cistron_is_mRNA]

		# Get boolean array of all cistrons that belong to at least one
		# polycistronic transcript
		polycistronic_cistron_indexes = []
		for rna_id in sim_data.process.transcription.rna_data['id']:
			cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes(rna_id)
			if len(cistron_indexes) > 1:
				polycistronic_cistron_indexes.extend(cistron_indexes)
		is_polycistronic = np.zeros(len(cistron_ids), bool)
		if len(polycistronic_cistron_indexes) > 0:
			is_polycistronic[np.array(list(set(polycistronic_cistron_indexes)))] = True
		mRNA_cistron_is_polycistronic = is_polycistronic[cistron_is_mRNA]

		# Get expression
		expected_mRNA_cistron_exp = sim_data.process.transcription.cistron_expression[CONDITION][cistron_is_mRNA]
		rna_exp = sim_data.process.transcription.rna_expression[CONDITION]
		actual_mRNA_cistron_exp = cistron_tu_mapping_matrix.dot(rna_exp)[cistron_is_mRNA]

		# Normalize expression
		expected_mRNA_cistron_exp /= expected_mRNA_cistron_exp.sum()
		actual_mRNA_cistron_exp /= actual_mRNA_cistron_exp.sum()

		# Find cistrons with more than a 10-fold difference between actual vs
		# expected expression
		exp_diff = np.log10(actual_mRNA_cistron_exp + NUMERICAL_ZERO) - np.log10(expected_mRNA_cistron_exp + NUMERICAL_ZERO)
		largest_diff_indexes = np.where(np.abs(exp_diff) > np.log10(LABEL_BOUNDARY))[0]

		plt.figure(figsize=(12, 12))
		ls = np.logspace(-8, -1, 8)
		plt.plot(ls, LABEL_BOUNDARY*ls, c='#dddddd', ls='--')
		plt.plot(ls, 1/LABEL_BOUNDARY*ls, c='#dddddd', ls='--')

		plt.scatter(
			expected_mRNA_cistron_exp, actual_mRNA_cistron_exp,
			c='#dddddd', s=2)

		# Highlight genes that are polycistronic
		plt.scatter(
			expected_mRNA_cistron_exp[mRNA_cistron_is_polycistronic],
			actual_mRNA_cistron_exp[mRNA_cistron_is_polycistronic],
			c='#333333', s=4, label='polycistronic')

		# Highlight cistrons that encode for RNAPs and ribosomal proteins
		plt.scatter(
			expected_mRNA_cistron_exp[mRNA_cistron_is_ribosomal_protein],
			actual_mRNA_cistron_exp[mRNA_cistron_is_ribosomal_protein],
			c='r', s=4, label='ribosomal proteins')
		plt.scatter(
			expected_mRNA_cistron_exp[mRNA_cistron_is_rnap],
			actual_mRNA_cistron_exp[mRNA_cistron_is_rnap],
			c='b', s=4, label='RNAP subunits')

		# Label cistrons that have more than a 10-fold difference
		for index in largest_diff_indexes:
			plt.text(
				expected_mRNA_cistron_exp[index],
				actual_mRNA_cistron_exp[index] * 1.05,
				mRNA_cistron_ids[index],
				ha='center', va='bottom', fontsize=5)

		plt.xlabel('Original mRNA cistron expression')
		plt.ylabel('Actual mRNA cistron expression')
		plt.xlim([1e-8, 1e-1])
		plt.ylim([1e-8, 1e-1])
		plt.xscale('log')
		plt.yscale('log')
		plt.legend()

		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
