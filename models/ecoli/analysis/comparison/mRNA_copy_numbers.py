"""
Generates a scatter plot of mRNA cistron copy numbers that are expected from the
expression levels calculated in the ParCa vs actual copies in the simulations
for a set of simulations run with operons, and a set run without operons.
"""
import itertools
import os
from typing import Tuple

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


FIGSIZE = (12, 6)
BOUNDS = [0, 2.5]
P_VALUE_THRESHOLD = 1e-3
N_BOOTSTRAP = 20000
NUMERICAL_ZERO = 1e-30

SEED = 0


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, validation_data2 = self.setup(input_sim_dir)

		np.random.seed(SEED)

		cell_paths = ap2.get_cells()
		simOutDir = os.path.join(cell_paths[0], 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Get mask for mRNA genes that are
		# i) Does not encode for ribosomal proteins or RNAPs
		# ii) not affected by manual overexpression/underexpression
		# to avoid confounding effects that come from manual adjustments of
		# expression levels.
		all_cistron_ids = sim_data2.process.transcription.cistron_data['id']
		is_mRNA = sim_data2.process.transcription.cistron_data['is_mRNA']
		assert np.all(
			mRNA_cistron_ids == all_cistron_ids[is_mRNA])

		mRNA_is_rnap_or_rprotein = np.logical_or(
				sim_data2.process.transcription.cistron_data['is_RNAP'],
				sim_data2.process.transcription.cistron_data['is_ribosomal_protein'])[is_mRNA]

		is_adjusted = np.zeros_like(is_mRNA, dtype=bool)
		all_rna_ids = sim_data2.process.transcription.rna_data['id']
		for adjusted_cistron_id in itertools.chain(
				sim_data2.adjustments.rna_expression_adjustments.keys(),
				sim_data2.adjustments.rna_deg_rates_adjustments.keys()):
			# Include cistrons whose expression is adjusted because they belong
			# to the same TU as the cistron that is bumped up
			adjusted_rna_indexes = sim_data2.process.transcription.cistron_id_to_rna_indexes(
				adjusted_cistron_id)
			adjusted_cistron_indexes = []
			for adjusted_rna_index in adjusted_rna_indexes:
				adjusted_cistron_indexes.extend(
					sim_data2.process.transcription.rna_id_to_cistron_indexes(
						all_rna_ids[adjusted_rna_index]))

			is_adjusted[adjusted_cistron_indexes] = True

		mRNA_is_adjusted = is_adjusted[is_mRNA]
		plot_mask = np.logical_and(~mRNA_is_rnap_or_rprotein, ~mRNA_is_adjusted)

		# Get mask for genes that are part of polycistronic transcription units
		polycistronic_cistron_indexes = []
		for rna_id in sim_data2.process.transcription.rna_data['id']:
			cistron_indexes = sim_data2.process.transcription.rna_id_to_cistron_indexes(rna_id)
			if len(cistron_indexes) > 1:
				polycistronic_cistron_indexes.extend(cistron_indexes)
		is_polycistronic = np.zeros(len(sim_data2.process.transcription.cistron_data), bool)
		if len(polycistronic_cistron_indexes) > 0:
			is_polycistronic[np.array(
				list(set(polycistronic_cistron_indexes)))] = True

		mRNA_mask_poly = is_polycistronic[is_mRNA]

		def read_data(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			# Sample initial mRNA counts from each cell
			all_initial_counts = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_cistron_counts',
				fun=lambda x: x[0])

			return all_initial_counts

		c1 = read_data(ap1)
		c2 = read_data(ap2)

		if len(c1) == 0 or len(c2) == 0:
			print('Skipping analysis -- not enough sims run.')
			return

		# Normalize counts from two conditions
		ratio = c1.mean(axis=0).sum()/c2.mean(axis=0).sum()
		c2 = ratio * c2.astype(np.float)

		# Calculate mean and variance
		m1 = c1.mean(axis=0)
		v1 = c1.var(axis=0)
		m2 = c2.mean(axis=0)
		v2 = c2.var(axis=0)

		def get_bootstrapped_samples(counts, sample_mean, combined_mean):
			# Bootstrap from counts corrected to have equal means
			n = counts.shape[0]

			bs_counts_mean = np.zeros((N_BOOTSTRAP, counts.shape[1]))
			bs_counts_var = np.zeros((N_BOOTSTRAP, counts.shape[1]))

			# Vectorizing this ran into memory issues with high values of
			# N_BOOTSTRAP
			for i in np.arange(N_BOOTSTRAP):
				bs_counts = (counts - sample_mean + combined_mean)[
					np.random.choice(n, n), :]

				bs_counts_mean[i, :] = bs_counts.mean(axis=0)
				bs_counts_var[i, :] = bs_counts.var(axis=0)

			return bs_counts_mean, bs_counts_var

		# Get bootstrapped samples from each count array
		n1 = c1.shape[0]
		n2 = c2.shape[1]
		combined_mean = (m1*n1 + m2*n2)/(n1 + n2)
		bm1, bv1 = get_bootstrapped_samples(c1, m1, combined_mean)
		bm2, bv2 = get_bootstrapped_samples(c2, m2, combined_mean)

		# Calculate t-score from the original sample and each of the
		# bootstrapped samples
		t_score = (m1 - m2)/np.sqrt(v1/n1 + v2/n2 + NUMERICAL_ZERO)
		t_score_bs = (bm1 - bm2)/np.sqrt(bv1/n1 + bv2/n2 + NUMERICAL_ZERO)

		# Calculate p-values of a two-tailed test using the empirical t-score
		# distribution calculated from the bootstrapped samples
		p_values = 2*np.minimum(
			((t_score <= t_score_bs).sum(axis=0) + 1)/(N_BOOTSTRAP + 1),
			((t_score >= t_score_bs).sum(axis=0) + 1)/(N_BOOTSTRAP + 1)
			)

		# Get mask for mRNAs with low p-values
		mRNA_mask_low_p = p_values < P_VALUE_THRESHOLD

		fig = plt.figure(figsize=FIGSIZE)

		for i, category in enumerate(['Polycistronic genes', 'Monocistronic genes']):
			ax = fig.add_subplot(1, 2, i + 1)
			ax.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
			category_mask = np.logical_xor(np.full_like(mRNA_mask_poly, i), mRNA_mask_poly)
			mask = np.logical_and.reduce((plot_mask, category_mask, ~mRNA_mask_low_p))
			ax.scatter(
				np.log10(m1[mask] + 1),
				np.log10(m2[mask] + 1),
				c='#cccccc', s=2, alpha=0.5,
				label=f'p ≥ {P_VALUE_THRESHOLD:g} (n = {mask.sum():d})',
				clip_on=False)
			# Highlight genes with low p-values
			mask = np.logical_and.reduce((plot_mask, category_mask, mRNA_mask_low_p))
			ax.scatter(
				np.log10(m1[mask] + 1),
				np.log10(m2[mask] + 1),
				c='r', s=2, alpha=0.5,
				label=f'p < {P_VALUE_THRESHOLD:g} (n = {mask.sum():d})',
				clip_on=False)

			ax.set_title(category)
			ax.set_xlabel('$\log_{10}$(mRNA copies + 1), old sims')
			ax.set_ylabel('$\log_{10}$(mRNA copies + 1), new sims')
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 20))
			ax.spines["left"].set_position(("outward", 20))
			ax.set_xlim(BOUNDS)
			ax.set_ylim(BOUNDS)
			ax.legend(loc=2)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Get bar plots of expression levels for operons with low p-values
		cistron_id_to_mRNA_index = {
			cistron_id: i for i, cistron_id in enumerate(mRNA_cistron_ids)
			}
		cistron_id_to_cistron_index = {
			cistron_id: i for i, cistron_id in enumerate(all_cistron_ids)
			}
		low_p_cistron_indexes = np.array([
			cistron_id_to_cistron_index[mRNA_cistron_ids[i]] for i
			in np.where(np.logical_and.reduce((plot_mask, mRNA_mask_poly, mRNA_mask_low_p)))[0]
			])

		operons_to_plot = [
			operon for operon in sim_data2.process.transcription.operons
			if np.any(np.isin(operon[0], low_p_cistron_indexes))]
		n_operons = len(operons_to_plot)

		# Get expression levels from each set
		operon_expression = []
		for operon_index, operon in enumerate(operons_to_plot):
			operon_cistron_ids = [all_cistron_ids[i] for i in operon[0]]
			operon_cistron_mRNA_indexes = np.array([
				cistron_id_to_mRNA_index[cistron_id]
				for cistron_id in operon_cistron_ids])

			operon_expression.append((
				m1[operon_cistron_mRNA_indexes], m2[operon_cistron_mRNA_indexes]
				))

		# Sort operons in descending order of average expression
		avg_expression = np.array([exp[0].mean() for exp in operon_expression])
		plot_order = np.argsort(avg_expression)[::-1]

		# Get mapping from cistron ids to gene names and transcription direction
		cistron_id_to_gene_name = {
			gene['cistron_id']: gene['symbol']
			for gene in sim_data2.process.replication.gene_data}
		cistron_is_forward = sim_data2.process.transcription.cistron_data['is_forward']

		fig = plt.figure()
		fig.set_size_inches(30, 4 * (n_operons//7 + 1))

		gs = gridspec.GridSpec(n_operons//7 + 1, 7)

		for i, operon_index in enumerate(plot_order):
			ax = plt.subplot(gs[i//7, i % 7])
			operon = operons_to_plot[operon_index]
			cistron_indexes_in_operon = operon[0]
			is_forward = cistron_is_forward[cistron_indexes_in_operon[0]]

			if not is_forward:
				cistron_indexes_in_operon = cistron_indexes_in_operon[::-1]

			# Cistron IDs with low p-values are starred
			operon_gene_names = [
				'*' + cistron_id_to_gene_name[all_cistron_ids[i]]
				if (i in low_p_cistron_indexes)
				else cistron_id_to_gene_name[all_cistron_ids[i]]
				for i in cistron_indexes_in_operon
				]
			operon_size = len(operon_gene_names)

			if is_forward:
				exp_without_operons = operon_expression[operon_index][0]
				exp_with_operons = operon_expression[operon_index][1]
			else:
				exp_without_operons = operon_expression[operon_index][0][::-1]
				exp_with_operons = operon_expression[operon_index][1][::-1]

			ax.bar(
				np.arange(operon_size) - 0.2,
				exp_without_operons,
				width=0.4, label='without operons')
			ax.bar(
				np.arange(operon_size) + 0.2,
				exp_with_operons,
				width=0.4, label='with operons')
			ax.set_xticks(np.arange(operon_size))
			ax.set_xticklabels(operon_gene_names, rotation=90)
			ax.set_ylabel('mRNA counts')

			if i == 0:
				ax.legend()

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_bar_plots', metadata)
		plt.close('all')


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == '__main__':
	Plot().cli()
