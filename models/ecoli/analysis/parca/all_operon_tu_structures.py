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
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id']
			for cistron in sim_data.process.transcription.cistron_data
			}
		cistron_coordinates = sim_data.process.transcription.cistron_data['replication_coordinate']
		tu_ids = [x[:-3] for x in sim_data.process.transcription.rna_data['id']]

		# Divide mapping matrix into individual operons
		cistron_indexes, tu_indexes = cistron_tu_mapping_matrix.nonzero()

		visited_cistron_indexes = set()
		visited_tu_indexes = set()
		all_operons = []

		def tu_DFS(tu_index, operon_cistron_indexes, operon_tu_indexes):
			"""
            Recursive function to look for indexes of TUs and cistrons that
            belong to the same operon as the TU with the given index.
            """
			visited_tu_indexes.add(tu_index)
			operon_tu_indexes.append(tu_index)

			for i in cistron_indexes[tu_indexes == tu_index]:
				if i not in visited_cistron_indexes:
					cistron_DFS(i, operon_cistron_indexes, operon_tu_indexes)

		def cistron_DFS(index, all_row_indexes, all_column_indexes):
			"""
            Recursive function to look for columns and rows in matrix A that should
            be grouped into the same NNLS problem.
            """
			visited_cistron_indexes.add(index)
			all_row_indexes.append(index)

			for i in tu_indexes[cistron_indexes == index]:
				if i not in visited_tu_indexes:
					tu_DFS(i, all_row_indexes, all_column_indexes)

		# Loop through each TU index
		for tu_index in range(len(tu_ids)):
			# Search for cistrons and TUs that can be grouped together into the
			# same operon
			if tu_index not in visited_tu_indexes:
				operon_cistron_indexes = []
				operon_tu_indexes = []
				tu_DFS(tu_index, operon_cistron_indexes, operon_tu_indexes)

				# Skip single-gene operons
				if len(operon_cistron_indexes) == 1 and len(operon_tu_indexes) == 1:
					continue

				# Sort cistron indexes by coordinates
				operon_cistron_indexes = sorted(
					operon_cistron_indexes,
					key = lambda i: cistron_coordinates[i])

				all_operons.append((
					operon_cistron_indexes, operon_tu_indexes
					))

		# Sort operons by size (number of cistrons)
		all_operons = sorted(
			all_operons,
			key = lambda operon: len(operon[0]), reverse=True)
		n_operons = len(all_operons)

		fig = plt.figure()
		fig.set_size_inches(50, 50)

		gs = gridspec.GridSpec(n_operons//30 + 1, 30)

		for (operon_index, (operon_cistron_indexes, operon_tu_indexes)) in enumerate(all_operons):
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
			ax.set_xticks(np.arange(len(operon_tu_indexes)))
			ax.set_yticks(np.arange(len(operon_cistron_indexes)))
			ax.set_xticklabels(
				[tu_ids[i] for i in operon_tu_indexes],
				rotation=90, fontsize=6)
			ax.set_yticklabels(
				[cistron_id_to_gene_id[cistron_ids[i]] for i in operon_cistron_indexes],
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
