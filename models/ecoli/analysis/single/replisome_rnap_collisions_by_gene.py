"""
Plots the number of collisions between RNAPs and replisomes that occur on each
gene. Only the top N genes with the most collisions are plotted.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2019
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import numpy as np
from six.moves import cPickle, range

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


PLOT_TOP_N_GENES = 25

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Load data
		n_headon = rnap_data_reader.readColumn("n_headon_collisions")
		n_codirectional = rnap_data_reader.readColumn("n_codirectional_collisions")
		headon_coordinates = rnap_data_reader.readColumn(
			"headon_collision_coordinates")
		codirectional_coordinates = rnap_data_reader.readColumn(
			"codirectional_collision_coordinates")

		# Get total number of collisions
		n_total_headon = n_headon.sum()
		n_total_codirectional = n_codirectional.sum()

		# Flatten collision coordinates
		headon_coordinates_flat = headon_coordinates[
			np.logical_not(np.isnan(headon_coordinates))].flatten()
		codirectional_coordinates_flat = codirectional_coordinates[
			np.logical_not(np.isnan(codirectional_coordinates))].flatten()

		# Get gene data from sim_data
		gene_ids = sim_data.process.transcription.rnaData["geneId"]
		gene_coordinates = sim_data.process.transcription.rnaData["replicationCoordinate"]
		gene_directions = sim_data.process.transcription.rnaData["direction"]
		gene_lengths = sim_data.process.transcription.rnaData["length"].asNumber(units.nt)

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Compute boundaries for each gene
		gene_boundaries = []
		for coord, direction, length in zip(
				gene_coordinates, gene_directions, gene_lengths):
			if direction:
				# Check for genes that loop around terC
				assert coord + length <= replichore_lengths[0]
				gene_boundaries.append((coord, coord + length))
			else:
				assert coord - length >= -replichore_lengths[1]
				gene_boundaries.append((coord - length, coord))

		# Initialize arrays
		n_codirectional_per_gene = np.zeros_like(gene_coordinates)
		n_headon_per_gene = np.zeros_like(gene_coordinates)

		# Identify which gene each collision occurred in
		for coord in codirectional_coordinates_flat:
			for i, boundaries in enumerate(gene_boundaries):
				if boundaries[0] < coord < boundaries[1]:
					n_codirectional_per_gene[i] += 1
					continue

		for coord in headon_coordinates_flat:
			for i, boundaries in enumerate(gene_boundaries):
				if boundaries[0] < coord < boundaries[1]:
					n_headon_per_gene[i] += 1
					continue

		# Sort by number of collisions
		codirectional_rank = np.argsort(n_codirectional_per_gene)[::-1][:PLOT_TOP_N_GENES]
		headon_rank = np.argsort(n_headon_per_gene)[::-1][:PLOT_TOP_N_GENES]

		# Get common names of top N genes
		codirectional_top_genes = [sim_data.common_names.genes[gene_ids[i]][0]
			for i in codirectional_rank]
		headon_top_genes = [sim_data.common_names.genes[gene_ids[i]][0]
			for i in headon_rank]

		# Plot
		plt.figure(figsize = (13, 3.5))

		ax = plt.subplot(1, 2, 1)
		ax.bar(list(range(PLOT_TOP_N_GENES)), n_codirectional_per_gene[codirectional_rank],
			color="darkblue")
		ax.set_xticks(list(range(PLOT_TOP_N_GENES)))
		ax.set_xticklabels(codirectional_top_genes, rotation=90)
		ax.set_title("Co-directional (Total = %d)"%(n_total_codirectional, ))
		ax.set_ylabel("Number of collisions")
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		ax = plt.subplot(1, 2, 2)
		ax.bar(list(range(PLOT_TOP_N_GENES)), n_headon_per_gene[headon_rank],
			color="crimson")
		ax.set_xticks(list(range(PLOT_TOP_N_GENES)))
		ax.set_xticklabels(headon_top_genes, rotation=90)
		ax.set_title("Head-on (Total = %d)"%(n_total_headon, ))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
