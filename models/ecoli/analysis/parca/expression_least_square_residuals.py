"""
Plots the differences bewteen the expected expression levels of each RNA cistron
vs the actual expression levels of each cistron calculated by multiplying the
cistron-TU mapping matrix with the expression levels of each transcription unit.
Any difference between the two are residuals from the nonnegative least squares
problem that was used to solve for the cistron expression levels.
"""

import pickle

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load data from sim_data
		all_conditions = list(sim_data.process.transcription.cistron_expression.keys())
		cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		cistron_is_ribosomal_protein = sim_data.process.transcription.cistron_data['is_ribosomal_protein']
		cistron_is_rnap = sim_data.process.transcription.cistron_data['is_RNAP']

		# Plot residuals for each condition
		fig = plt.figure()
		n_plots = len(all_conditions)
		fig.set_size_inches(23, 5.5 * np.ceil(n_plots/4))
		gs = gridspec.GridSpec(int(np.ceil(n_plots/4)), 4)

		for (i, condition) in enumerate(all_conditions):
			condition_plt = plt.subplot(gs[i//4, i%4])
			expected_cistron_exp = sim_data.process.transcription.cistron_expression[condition]
			rna_exp = sim_data.process.transcription.rna_expression[condition]
			actual_cistron_exp = cistron_tu_mapping_matrix.dot(rna_exp)

			condition_plt.scatter(
				expected_cistron_exp, actual_cistron_exp, c='#bbbbbb', s=2)

			# Highlight genes that encode for RNAPs and ribosomal proteins
			condition_plt.scatter(
				expected_cistron_exp[cistron_is_ribosomal_protein],
				actual_cistron_exp[cistron_is_ribosomal_protein], c='r', s=3,
				label='ribosomal proteins')
			condition_plt.scatter(
				expected_cistron_exp[cistron_is_rnap],
				actual_cistron_exp[cistron_is_rnap], c='b', s=3,
				label='RNAP subunits')

			condition_plt.set_title(condition)
			condition_plt.set_xlabel('Expected cistron expression (b)')
			condition_plt.set_ylabel('Actual cistron expression (Ax)')
			condition_plt.set_xscale('log')
			condition_plt.set_yscale('log')
			if i == 0:
				condition_plt.legend()

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
