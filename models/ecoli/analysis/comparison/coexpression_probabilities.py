"""
Compare coexpression probabilities of genes in the same operon at the mRNA and
protein level.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from scipy import stats

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError


FIGSIZE = (14.85, 4)
LOW_EXP_THRESHOLD = 100

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		# Check two sims have same list of genes
		if np.any(sim_data1.process.transcription.cistron_data['id'] != sim_data2.process.transcription.cistron_data['id']):
			print('Skipping analysis -- two sims must have same set of genes.')
			return

		# Get list of gene groups that are part of the same operon
		is_mRNA = sim_data2.process.transcription.cistron_data['is_mRNA']
		gene_group_indexes = [
			operon[0] for operon in sim_data2.process.transcription.operons
			if len(operon[0]) > 1 and np.all(is_mRNA[operon[0]])]
		gene_index_to_mRNA_index = {
			gene_index: mRNA_index
			for (mRNA_index, gene_index) in enumerate(np.where(is_mRNA)[0])
			}
		mRNA_group_indexes = [
			[gene_index_to_mRNA_index[x] for x in gene_indexes]
			for gene_indexes in gene_group_indexes
			]
		cistron_id_to_cistron_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(sim_data2.process.transcription.cistron_data['id'])
			}
		cistron_index_to_monomer_index = {
			cistron_id_to_cistron_index[cistron_id]: monomer_index
			for (monomer_index, cistron_id)
			in enumerate(sim_data2.process.translation.monomer_data['cistron_id'])
			}
		monomer_group_indexes = [
			[cistron_index_to_monomer_index[x] for x in gene_indexes]
			for gene_indexes in gene_group_indexes
			]

		def read_sims_for_mRNA(ap):
			mRNA_coexp_probs = np.zeros(len(gene_group_indexes))
			is_expressed = np.zeros(len(gene_group_indexes), dtype=np.bool)

			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			all_mRNA_counts = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_cistron_counts', ignore_exception=True)

			for i in range(len(gene_group_indexes)):
				operon_mRNA_counts = all_mRNA_counts[:, np.array(mRNA_group_indexes[i])]
				mRNA_coexp_probs[i] = (
					(np.all(operon_mRNA_counts, axis=1).sum() + 1) /
					(np.any(operon_mRNA_counts, axis=1).sum() + 1)
				)
				is_expressed[i] = np.all(operon_mRNA_counts.sum(axis=0))

			return mRNA_coexp_probs, is_expressed

		def read_sims_for_protein(ap):
			protein_coexp_probs = np.zeros(len(gene_group_indexes))
			is_lowly_expressed = np.zeros(len(gene_group_indexes), dtype=np.bool)

			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			all_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)

			for i in range(len(gene_group_indexes)):

				operon_monomer_counts = all_monomer_counts[:, np.array(monomer_group_indexes[i])]
				protein_coexp_probs[i] = (
					(np.all(operon_monomer_counts, axis=1).sum() + 1) /
					(np.any(operon_monomer_counts, axis=1).sum() + 1)
				)
				is_lowly_expressed[i] = (operon_monomer_counts.mean() < LOW_EXP_THRESHOLD)

			return protein_coexp_probs, is_lowly_expressed


		mRNA_p1, is_expressed = read_sims_for_mRNA(ap1)
		mRNA_p2, _ = read_sims_for_mRNA(ap2)
		protein_p1, is_lowly_expressed = read_sims_for_protein(ap1)
		protein_p2, _ = read_sims_for_protein(ap2)

		fig = plt.figure(figsize=FIGSIZE)
		gs = fig.add_gridspec(
			2, 6, width_ratios=(4, 1, 4, 1, 4, 1), height_ratios=(1, 4))

		def draw_plot(p1, p2, grid_i, grid_j, y_max):
			scatter_ax = fig.add_subplot(gs[grid_i, grid_j])
			scatter_ax.plot([0, 1], [0, 1], ls='--', lw=2, c='k', alpha=0.1)
			scatter_ax.scatter(
				p1, p2,
				alpha=0.5, s=5, c='k', clip_on=False, edgecolors='none')

			scatter_ax.set_xlim([0, 1])
			scatter_ax.set_ylim([0, 1])
			scatter_ax.set_xlabel('Reference')
			scatter_ax.set_ylabel('Input')

			scatter_ax.spines["top"].set_visible(False)
			scatter_ax.spines["right"].set_visible(False)
			scatter_ax.spines["bottom"].set_position(("outward", 20))
			scatter_ax.spines["left"].set_position(("outward", 20))

			x = np.linspace(0, 1, 1000)
			kde1 = stats.gaussian_kde(p1)
			kde2 = stats.gaussian_kde(p2)

			hist1_ax = fig.add_subplot(gs[grid_i - 1, grid_j], sharex=scatter_ax)
			hist1_ax.fill_between(x, kde1(x), alpha=0.5)
			hist1_ax.set_xlim([0, 1])
			hist1_ax.set_ylim([0, y_max])
			hist1_ax.set_yticks([])
			hist1_ax.spines["top"].set_visible(False)
			hist1_ax.spines["right"].set_visible(False)
			hist1_ax.spines["left"].set_visible(False)
			hist1_ax.spines["bottom"].set_visible(False)
			plt.setp(hist1_ax.get_xaxis(), visible=False)

			hist2_ax = fig.add_subplot(gs[grid_i, grid_j + 1], sharey=scatter_ax)
			hist2_ax.fill_betweenx(x, kde2(x), fc='C1', alpha=0.5)
			hist2_ax.set_ylim([0, 1])
			hist2_ax.set_xlim([0, y_max])
			hist2_ax.set_xticks([])
			hist2_ax.spines["top"].set_visible(False)
			hist2_ax.spines["right"].set_visible(False)
			hist2_ax.spines["left"].set_visible(False)
			hist2_ax.spines["bottom"].set_visible(False)
			plt.setp(hist2_ax.get_yaxis(), visible=False)

		draw_plot(
			mRNA_p1[is_expressed],
			mRNA_p2[is_expressed],
			1, 0, 5)
		draw_plot(
			protein_p1[np.logical_and(is_expressed, is_lowly_expressed)],
			protein_p2[np.logical_and(is_expressed, is_lowly_expressed)],
			1, 2, 4)
		draw_plot(
			protein_p1[np.logical_and(is_expressed, ~is_lowly_expressed)],
			protein_p2[np.logical_and(is_expressed, ~is_lowly_expressed)],
			1, 4, 12)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
