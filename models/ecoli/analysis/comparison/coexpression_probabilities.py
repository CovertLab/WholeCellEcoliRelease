"""
Compare coexpression probabilities of operons.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError


FIGSIZE = (8, 4)

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

		def read_sims(ap):
			mRNA_coexp_probs = np.zeros(len(gene_group_indexes))
			protein_coexp_probs = np.zeros(len(gene_group_indexes))
			all_genes_expressed = np.zeros(len(gene_group_indexes), dtype=np.bool)

			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			all_mRNA_counts = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_cistron_counts')
			all_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts')
			n_timesteps = all_mRNA_counts.shape[0]

			for i in range(len(gene_group_indexes)):
				operon_mRNA_counts = all_mRNA_counts[:, np.array(mRNA_group_indexes[i])]
				mRNA_coexp_probs[i] = np.all(operon_mRNA_counts, axis=1).sum()/n_timesteps
				operon_monomer_counts = all_monomer_counts[:, np.array(monomer_group_indexes[i])]
				protein_coexp_probs[i] = np.all(operon_monomer_counts, axis=1).sum() / n_timesteps
				all_genes_expressed[i] = np.all(operon_mRNA_counts.sum(axis=0))

			return mRNA_coexp_probs, protein_coexp_probs, all_genes_expressed

		mRNA_probs1, protein_probs1, all_genes_expressed = read_sims(ap1)
		mRNA_probs2, protein_probs2, _ = read_sims(ap2)

		fig = plt.figure(figsize=FIGSIZE)
		ax0 = fig.add_subplot(1, 2, 1)
		ax0.plot([0, 1], [0, 1], ls='--', lw=2, c='k', alpha=0.05)
		ax0.scatter(
			mRNA_probs1[all_genes_expressed], mRNA_probs2[all_genes_expressed],
			alpha=0.5, s=5, c='k', clip_on=False, edgecolors='none')

		ax0.set_xlim([0, 1])
		ax0.set_ylim([0, 1])
		ax0.set_title('mRNA')
		ax0.set_xlabel('Reference')
		ax0.set_ylabel('Input')

		ax0.spines["top"].set_visible(False)
		ax0.spines["right"].set_visible(False)
		ax0.spines["bottom"].set_position(("outward", 20))
		ax0.spines["left"].set_position(("outward", 20))

		ax1 = fig.add_subplot(1, 2, 2)
		ax1.plot([0, 1], [0, 1], ls='--', lw=2, c='k', alpha=0.05)
		ax1.scatter(
			protein_probs1[all_genes_expressed], protein_probs2[all_genes_expressed],
			alpha=0.5, s=5, c='k', clip_on=False, edgecolors='none')

		ax1.set_xlim([0, 1])
		ax1.set_ylim([0, 1])
		ax1.set_title('protein')
		ax1.set_xlabel('Reference')
		ax1.set_ylabel('Input')

		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 20))
		ax1.spines["left"].set_position(("outward", 20))

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		import ipdb; ipdb.set_trace()

	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
