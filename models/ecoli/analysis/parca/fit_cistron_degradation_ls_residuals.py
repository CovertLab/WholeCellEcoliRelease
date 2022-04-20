"""
Plots the differences bewteen the expected degradation rates of each RNA cistron
vs the actual degradation rates of each cistron calculated from the fit
degradation rates of RNAs under basal minimum media. Any differences between the
two are residuals from the nonnegative least squares problem that was used to
solve for RNA-level degradation rates in sim_data.
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load data from sim_data
		rna_deg_rates = sim_data.process.transcription.rna_data['deg_rate'].asNumber(1/units.s)
		cistron_deg_rates = sim_data.process.transcription.cistron_data['deg_rate'].asNumber(1/units.s)
		rna_expression = sim_data.process.transcription.rna_expression['basal']
		all_cistron_ids = sim_data.process.transcription.cistron_data['id']

		# Build the relative abundance matrix between transcription units and
		# constituent cistrons
		cistron_indexes = []
		rna_indexes = []
		v = []

		for cistron_index, cistron_id in enumerate(all_cistron_ids):
			rna_indexes_this_cistron = sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id)
			v_this_cistron = np.zeros(len(rna_indexes_this_cistron))
			for i, rna_index in enumerate(rna_indexes_this_cistron):
				cistron_indexes.append(cistron_index)
				rna_indexes.append(rna_index)
				v_this_cistron[i] = rna_expression[rna_index]

			if v_this_cistron.sum() == 0:
				v_this_cistron[:] = 1./len(v_this_cistron)
			else:
				v_this_cistron = v_this_cistron/v_this_cistron.sum()

			v.extend(v_this_cistron)

		cistron_indexes = np.array(cistron_indexes)
		rna_indexes = np.array(rna_indexes)
		v = np.array(v)
		shape = (cistron_indexes.max() + 1, rna_indexes.max() + 1)

		cistron_tu_relative_abundancy_matrix = csr_matrix(
			(v, (cistron_indexes, rna_indexes)),
			shape=shape)

		# Calculate effective degradation rates of each cistron using relative
		# abundancy matrix
		effective_deg_rates = cistron_tu_relative_abundancy_matrix.dot(rna_deg_rates)

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

		# Get mask for cistrons that have one or more encoding TUs with
		# measured degradation rates
		measured_deg_rate_mask = (sim_data.process.transcription.cistron_tu_mapping_matrix.dot(
			sim_data.process.transcription.rna_data['deg_rate_is_measured'])).astype(np.bool)

		# Plot expected vs actual deg rates
		plt.figure(figsize=(8, 8))
		plt.scatter(
			cistron_deg_rates[np.logical_and(is_polycistronic, measured_deg_rate_mask)],
			effective_deg_rates[np.logical_and(is_polycistronic, measured_deg_rate_mask)],
			c='#bbbbbb', s=5, label='TU-level rates used')
		plt.scatter(
			cistron_deg_rates[np.logical_and(is_polycistronic, ~measured_deg_rate_mask)],
			effective_deg_rates[np.logical_and(is_polycistronic, ~measured_deg_rate_mask)],
			c='k', s=5, label='Inferred from gene-level rates')
		plt.xlabel('Expected k_deg (1/s)')
		plt.ylabel('Actual effective k_deg (1/s)')
		plt.xlim([2e-4, 0.04])
		plt.ylim([2e-4, 0.04])
		plt.xscale('log')
		plt.yscale('log')
		plt.legend()
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
