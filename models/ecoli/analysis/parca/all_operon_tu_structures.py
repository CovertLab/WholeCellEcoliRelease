"""
Plots the transcription unit structures of all operons being simulated in the
model as binary heatmaps. The cistrons within each operon are sorted by their
genomic coordinates. Monocistronic operons are not shown.
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
		cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		cistron_ids = sim_data.process.transcription.cistron_data['id']
		cistron_id_to_gene_name = {
			gene['cistron_id']: gene['symbol']
			for gene in sim_data.process.replication.gene_data
			}

		# Divide mapping matrix into individual operons
		cistron_indexes, tu_indexes = cistron_tu_mapping_matrix.nonzero()

		# Get polycistronic operons
		operons = [
			operon for operon in sim_data.process.transcription.operons
			if len(operon[0]) > 1
			]

		# Sort operons by size (number of cistrons)
		operons_sorted = sorted(
			operons, key = lambda operon: len(operon[0]), reverse=True)
		n_operons = len(operons_sorted)

		fig = plt.figure()
		fig.set_size_inches(50, 50)

		gs = gridspec.GridSpec(n_operons//30 + 1, 30)

		for (operon_index, (operon_cistron_indexes, operon_tu_indexes)) in enumerate(operons_sorted):
			ax = plt.subplot(gs[operon_index // 30, operon_index % 30])

			# Get TU structure array of this operon
			tu_structure = np.full(
				(len(operon_cistron_indexes), len(operon_tu_indexes)), 255)
			mask = np.isin(cistron_indexes, operon_cistron_indexes)
			for (i, j, v) in zip(
					cistron_indexes[mask],
					tu_indexes[mask],
					cistron_tu_mapping_matrix.data[mask]):
				tu_structure[
					np.where(operon_cistron_indexes == i)[0][0],
					np.where(operon_tu_indexes == j)[0][0]] = 0

			# Display TU structure as binary heatmap
			ax.imshow(tu_structure, cmap='gray')
			ax.set_xticks([])
			ax.set_yticks(np.arange(len(operon_cistron_indexes)))
			ax.set_yticklabels(
				[cistron_id_to_gene_name[cistron_ids[i]] for i in operon_cistron_indexes],
				fontsize=6)

			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['left'].set_visible(False)
			ax.tick_params(axis='both', bottom=False, left=False, pad=1)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
