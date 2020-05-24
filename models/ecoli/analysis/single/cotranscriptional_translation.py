'''
Patterns of co-transcriptional translation in the E. coli WCM

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/17/20
'''

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath, units

PLOT_TOP_N_GENES = 30


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception(
				'simOutDir does not currently exist as a directory')

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')
		all_mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		full_mRNA_counts = mRNA_counts_reader.readColumn('full_mRNA_counts')
		partial_mRNA_counts = mRNA_counts_reader.readColumn('partial_mRNA_counts')
		protein_ids = ribosome_reader.readAttribute('monomerIds')
		all_ribosome_counts = ribosome_reader.readColumn('n_ribosomes_per_transcript')
		chromosome_bound_ribosome_counts = ribosome_reader.readColumn(
			'n_ribosomes_on_partial_mRNA_per_transcript')

		# Load mappings from mRNA, protein ids to gene ids
		mRNA_id_to_gene_id = {
			mRNA_id: gene_id for (mRNA_id, gene_id)
			in zip(sim_data.process.transcription.rnaData['id'],
				sim_data.process.transcription.rnaData['geneId'])
			}
		protein_id_to_gene_id = {
			protein_id: mRNA_id_to_gene_id[mRNA_id] for (protein_id, mRNA_id)
			in zip(sim_data.process.translation.monomerData['id'],
				sim_data.process.translation.monomerData['rnaId'])
			}

		fig = plt.figure()
		fig.set_size_inches(8, 28)

		gs = gridspec.GridSpec(8, 1)

		# Plot counts of full/partial transcripts over time
		ax = plt.subplot(gs[0, 0])
		ax.plot(time, all_mRNA_counts.sum(axis=1), label='All mRNAs')
		ax.plot(time, full_mRNA_counts.sum(axis=1), label='Fully transcribed mRNAs')
		ax.plot(time, partial_mRNA_counts.sum(axis=1), label='Partially transcribed mRNAs')
		ax.legend()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_xlabel('Time [s]')
		ax.set_ylabel('Total Counts')

		# Get top N genes with the highest average count of active partially
		# transcribed mRNAs
		partial_mRNA_counts_mean = partial_mRNA_counts.mean(axis=0)
		rank = np.argsort(partial_mRNA_counts_mean)[::-1][:PLOT_TOP_N_GENES]
		ranked_genes = [sim_data.common_names.genes[
			mRNA_id_to_gene_id[mRNA_ids[i]]][0] for i in rank]

		ax = plt.subplot(gs[1, 0])
		ax.bar(range(PLOT_TOP_N_GENES), partial_mRNA_counts_mean[rank])
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(ranked_genes, rotation=90)
		ax.set_title(
			'Genes with highest average counts of partially transcribed active transcripts')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylabel('Average counts of\npartially transcribed transcripts')

		# Get top N genes with the highest proportions of active mRNAs that are
		# partially transcribed (minimum 10 average total transcripts)
		all_mRNA_counts_mean = all_mRNA_counts.mean(axis=0)
		partial_mRNA_proportions_all_mRNAs = partial_mRNA_counts_mean.sum()/all_mRNA_counts_mean.sum()
		partial_mRNA_proportions = np.nan_to_num(
			partial_mRNA_counts_mean/all_mRNA_counts_mean)
		partial_mRNA_proportions[all_mRNA_counts_mean < 10] = 0
		rank = np.argsort(partial_mRNA_proportions)[::-1][:PLOT_TOP_N_GENES]
		ranked_genes = [sim_data.common_names.genes[
			mRNA_id_to_gene_id[mRNA_ids[i]]][0] for i in rank]

		ax = plt.subplot(gs[2, 0])
		ax.bar(range(PLOT_TOP_N_GENES), partial_mRNA_proportions[rank])
		ax.axhline(
			partial_mRNA_proportions_all_mRNAs, ls='--', lw=2, color='k',
			label="All mRNAs: %.2f" % (partial_mRNA_proportions_all_mRNAs, )
			)
		ax.legend()
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(ranked_genes, rotation=90)
		ax.set_title('Genes with highest proportions of partially transcribed active transcripts\n($\geq 10$ total transcripts)')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylim([0, 1])
		ax.set_ylabel('Percentage of\npartially transcribed transcripts')

		# Get scatter plot of observed mRNA proportions vs the proportions
		# expected from the degradation rates and lengths of each transcript.
		# Assuming a constant transcription rate and a first-order degradation
		# rate, the distribution of ages of mRNA molecules follows an
		# exponential distribution, similar to CSTR reactors. Thus, the
		# proportions of partially transcribed transcripts can be estimated by
		# taking an integral of this distribution function from 0 to T, where T
		# is the length of time required to produce a full transcript.
		is_mRNA = sim_data.process.transcription.rnaData['isMRna']
		mRNA_lengths = sim_data.process.transcription.rnaData['length'][is_mRNA].asNumber(units.nt)
		mRNA_deg_rates = sim_data.process.transcription.rnaData['degRate'][is_mRNA].asNumber(1/units.s)
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		elongation_rate = sim_data.process.transcription.rnaPolymeraseElongationRateDict[
			nutrients].asNumber(units.nt/units.s)

		mRNA_transcription_time = mRNA_lengths.astype(np.float64)/elongation_rate
		expected_partial_mRNA_proportions = np.ones_like(mRNA_deg_rates) - np.exp(
			-np.multiply(mRNA_deg_rates, mRNA_transcription_time))
		mask = all_mRNA_counts_mean > 10  # Plot only the highly-expressed mRNAs

		ax = plt.subplot(gs[3:5, 0])
		ax.scatter(
			expected_partial_mRNA_proportions[mask],
			partial_mRNA_proportions[mask],
			marker='o', s=20, alpha=0.3, color='k'
			)
		ax.plot([0, 1], [0, 1], '--r')
		ax.set_title(
			'Validation of partially transcribed transcript proportions with expected values')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_xlim([0, 1])
		ax.set_ylim([0, 1])
		ax.set_xlabel('Expected proportions of partially transcribed transcripts')
		ax.set_ylabel('Observed proportions of partially transcribed transcripts')

		# Plot counts of all active ribosomes/ribosomes that are translating
		# partially transcribed mRNAs over time
		ax = plt.subplot(gs[5, 0])
		free_ribosome_counts = all_ribosome_counts - chromosome_bound_ribosome_counts
		ax.plot(time, all_ribosome_counts.sum(axis=1), label='All ribosomes')
		ax.plot(time, free_ribosome_counts.sum(axis=1), label='Free ribosomes')
		ax.plot(time, chromosome_bound_ribosome_counts.sum(axis=1), label='Chromosome-bound ribosomes')
		ax.legend()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_xlabel('Time [s]')
		ax.set_ylabel('Total Counts')

		# Get top N genes with the highest average counts active ribosomes that
		# are bound to the chromosome
		bound_ribosome_counts_mean = chromosome_bound_ribosome_counts.mean(axis=0)
		rank = np.argsort(bound_ribosome_counts_mean)[::-1][:PLOT_TOP_N_GENES]
		ranked_genes = [sim_data.common_names.genes[
			protein_id_to_gene_id[protein_ids[i]]][0] for i in rank]

		ax = plt.subplot(gs[6, 0])
		ax.bar(range(PLOT_TOP_N_GENES), bound_ribosome_counts_mean[rank])
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(ranked_genes, rotation=90)
		ax.set_title(
			'Genes with highest average counts of chromosome-bound ribosomes')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylabel('Average counts of\nchromosome-bound ribosomes')

		# Get top N genes with the highest proportions of active ribosomes that
		# are bound to the chromosome (minimum 10 average total ribosomes)
		all_ribosome_counts_mean = all_ribosome_counts.mean(axis=0)
		bound_ribosome_counts_mean = chromosome_bound_ribosome_counts.mean(axis=0)
		bound_ribosome_proportions_all_ribosomes = bound_ribosome_counts_mean.sum() / all_ribosome_counts_mean.sum()
		bound_ribosome_proportions = np.nan_to_num(
			bound_ribosome_counts_mean / all_ribosome_counts_mean)
		bound_ribosome_proportions[all_ribosome_counts_mean < 10] = 0
		rank = np.argsort(bound_ribosome_proportions)[::-1][:PLOT_TOP_N_GENES]
		ranked_genes = [sim_data.common_names.genes[
			protein_id_to_gene_id[protein_ids[i]]][0] for i in rank]

		ax = plt.subplot(gs[7, 0])
		ax.bar(range(PLOT_TOP_N_GENES), bound_ribosome_proportions[rank])
		ax.axhline(
			bound_ribosome_proportions_all_ribosomes, ls='--', lw=2, color='k',
			label="All ribosomes: %.2f" % (bound_ribosome_proportions_all_ribosomes, )
			)
		ax.legend()
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(ranked_genes, rotation=90)
		ax.set_title(
			'Genes with highest proportions of chromosome-bound ribosomes\n($\geq 10$ total ribosomes)')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylim([0, 1])
		ax.set_ylabel('Proportions of\nchromosome-bound ribosomes')

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
