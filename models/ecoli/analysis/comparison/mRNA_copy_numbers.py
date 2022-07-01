"""
Generates a comparison scatter plot of mRNA copy numbers from two sets of
simulations, one without operons and one with operons.
"""
import itertools
import os
from typing import Tuple

import csv
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


GENERATE_OPERON_TABLE = False
FIGSIZE = (8, 4.1)
BOUNDS = [[0, 2], [0, 2.5]]
SELECTION_RATIO = 0.1
NUMERICAL_ZERO = 1e-30
OPERON_COUNT_CUTOFF = 10

EVIDENCE_CODE_TO_DESCRIPTIONS = {
	'EV-COMP-HINF': 'Human inference from comp. ev.',
	'EV-COMP-AINF': 'Automated inference from comp. ev.',
	'EV-EXP-IEP': 'Expression pattern',
	'EV-EXP-IEP-COREGULATION': 'Coregulation',
	'EV-IC-ADJ-GENES-SAME-BIO-PROCESS': 'Gene adjacency',
	'EV-IC': 'Curator inference',
	'EV-EXP-IDA-TRANSCRIPT-LEN-DETERMINATION': 'Transcript lengths',
	'EV-EXP-IDA-BOUNDARIES-DEFINED': 'Transcription boundaries',
	'EV-EXP-IMP-POLAR-MUTATION': 'Polar effects',
	0: 'No evidence',
	1: 'Single evidence',
	2: 'Multiple evidences',
	}


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

		n1 = c1.shape[0]
		n2 = c2.shape[0]

		# Calculate t-score and get cutoff
		abs_t_scores = np.abs(m1 - m2)/(np.sqrt(v1/n1 + v2/n2) + NUMERICAL_ZERO)
		n_plotted_poly_genes = np.logical_and(plot_mask, mRNA_mask_poly).sum()
		t_score_cutoff = np.sort(
			abs_t_scores[np.logical_and(plot_mask, mRNA_mask_poly)]
			)[-int(n_plotted_poly_genes * SELECTION_RATIO)]

		# Get mask for mRNAs with high absolute t-scores
		mRNA_mask_high_t = (abs_t_scores >= t_score_cutoff)

		fig = plt.figure(figsize=FIGSIZE)

		for i, category in enumerate(['Polycistronic genes', 'Monocistronic genes']):
			ax = fig.add_subplot(1, 2, i + 1)
			ax.plot(BOUNDS[i], BOUNDS[i], ls='--', lw=2, c='k', alpha=0.05)
			category_mask = np.logical_xor(np.full_like(mRNA_mask_poly, i), mRNA_mask_poly)
			mask = np.logical_and.reduce((plot_mask, category_mask, ~mRNA_mask_high_t))
			ax.scatter(
				np.log10(m1[mask] + 1),
				np.log10(m2[mask] + 1),
				c='#555555', edgecolor='none', s=10, alpha=0.4,
				label=f'|t| < {t_score_cutoff:.1f} (n = {mask.sum():d})',
				clip_on=False)
			# Highlight genes with high t-scores
			mask = np.logical_and.reduce((plot_mask, category_mask, mRNA_mask_high_t))
			ax.scatter(
				np.log10(m1[mask] + 1),
				np.log10(m2[mask] + 1),
				c='C3', edgecolor='none', s=10, alpha=0.4,
				label=f'|t| ≥ {t_score_cutoff:.1f} (n = {mask.sum():d})',
				clip_on=False)

			ax.set_title(category)
			ax.set_xlabel('$\log_{10}$(mRNA copies + 1), old sims')
			ax.set_ylabel('$\log_{10}$(mRNA copies + 1), new sims')
			ax.set_xticks(np.arange(BOUNDS[i][0], BOUNDS[i][1] + 0.5, 0.5))
			ax.set_yticks(np.arange(BOUNDS[i][0], BOUNDS[i][1] + 0.5, 0.5))
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim(BOUNDS[i])
			ax.set_ylim(BOUNDS[i])
			ax.legend(loc=2, prop={'size': 8})

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Get bar plots of expression levels for operons with high t-scores
		cistron_id_to_mRNA_index = {
			cistron_id: i for i, cistron_id in enumerate(mRNA_cistron_ids)
			}
		cistron_id_to_cistron_index = {
			cistron_id: i for i, cistron_id in enumerate(all_cistron_ids)
			}
		high_t_cistron_indexes = np.array([
			cistron_id_to_cistron_index[mRNA_cistron_ids[i]] for i
			in np.where(np.logical_and.reduce((plot_mask, mRNA_mask_poly, mRNA_mask_high_t)))[0]
			])

		high_t_operons = [
			operon for operon in sim_data2.process.transcription.operons
			if np.any(np.isin(operon[0], high_t_cistron_indexes))]
		n_high_t_operons = len(high_t_operons)

		# Get expression levels from each set
		operon_expression = []
		for operon_index, operon in enumerate(high_t_operons):
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
		fig.set_size_inches(30, 4 * (n_high_t_operons//7 + 1))

		gs = gridspec.GridSpec(n_high_t_operons//7 + 1, 7)

		for i, operon_index in enumerate(plot_order):
			ax = plt.subplot(gs[i//7, i % 7])
			operon = high_t_operons[operon_index]
			cistron_indexes_in_operon = operon[0]
			is_forward = cistron_is_forward[cistron_indexes_in_operon[0]]

			if not is_forward:
				cistron_indexes_in_operon = cistron_indexes_in_operon[::-1]

			# Cistron IDs with high t-scores are starred
			operon_gene_names = [
				'*' + cistron_id_to_gene_name[all_cistron_ids[i]]
				if (i in high_t_cistron_indexes)
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

		# (Optional) Get table of operons with high t-score genes
		if GENERATE_OPERON_TABLE:
			with open(os.path.join(plotOutDir, plotOutFileName + '_operon_table.tsv'), 'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow([
					'operon_name', 'first_gene', 'last_gene', 'max_t'
					])

				for operon_index in plot_order:
					operon = high_t_operons[operon_index]
					cistron_indexes = operon[0]
					t_scores_this_operon = [
						abs_t_scores[cistron_id_to_mRNA_index[all_cistron_ids[i]]]
						for i in cistron_indexes]

					is_forward = cistron_is_forward[cistron_indexes[0]]
					if is_forward:
						gene_list = [cistron_id_to_gene_name[all_cistron_ids[i]] for i in cistron_indexes]
					else:
						gene_list = [cistron_id_to_gene_name[all_cistron_ids[i]] for i in cistron_indexes[::-1]]

					writer.writerow([
						'-'.join(gene_list),
						cistron_id_to_gene_name[all_cistron_ids[cistron_indexes[0]]],
						cistron_id_to_gene_name[all_cistron_ids[cistron_indexes[-1]]],
						max(t_scores_this_operon),
						])

		# Get violin plots of log t-scores of each evidence code
		# Map each operon to the list of evidence codes for the transcription
		# units in the operon
		all_operons = sim_data2.process.transcription.operons
		corrected_cistron_indexes = np.where(
			sim_data2.process.transcription.cistron_data['uses_corrected_seq_counts']
			)[0]
		operon_index_to_evidence_codes = {}

		for (i, operon) in enumerate(all_operons):
			# Skip monocistronic operons
			if len(operon[0]) == 1:
				continue

			# Skip operons that contain cistrons whose expression was corrected
			if np.any(np.isin(operon[0], corrected_cistron_indexes)):
				continue

			evidence_codes = []
			for rna_index in operon[1]:
				evidence_codes.extend(
					sim_data2.process.transcription.rna_id_to_evidence_codes.get(
						all_rna_ids[rna_index][:-3], [])
					)

			operon_index_to_evidence_codes[i] = list(set(evidence_codes))

		# Map evidence code and multiplicity to maximum log t-scores of each
		# operon
		all_operons_log_t_scores = []
		evidence_code_to_log_t_scores = {}
		evidence_multiplicity_to_log_t_scores = {}

		for (operon_index, evidence_codes) in operon_index_to_evidence_codes.items():
			operon = all_operons[operon_index]
			operon_cistron_ids = [all_cistron_ids[i] for i in operon[0]]
			operon_cistron_mRNA_indexes = np.array([
				cistron_id_to_mRNA_index[cistron_id]
				for cistron_id in operon_cistron_ids])
			max_log_t = np.log10(
				abs_t_scores[operon_cistron_mRNA_indexes].max() + 1)

			for evidence_code in evidence_codes:
				evidence_code_to_log_t_scores.setdefault(
					evidence_code, []).append(max_log_t)

			# Evidence multiplicity is maxed out at two
			evidence_multiplicity_to_log_t_scores.setdefault(
				min(len(evidence_codes), 2), []).append(max_log_t)

			all_operons_log_t_scores.append(max_log_t)

		def sort_t_score_dict(t_score_dict):
			sorted_t_score_dict = {}

			for key, value in t_score_dict.items():
				if len(value) < OPERON_COUNT_CUTOFF:
					continue
				sorted_t_score_dict[key] = np.array(value)

			sorted_t_score_dict = dict(sorted(
				sorted_t_score_dict.items(), key=lambda item: np.median(item[1]),
				reverse=True))

			return sorted_t_score_dict

		all_operons_log_t_scores = {
			'All operons': np.array(all_operons_log_t_scores)
			}
		evidence_code_to_log_t_scores = sort_t_score_dict(
			evidence_code_to_log_t_scores)
		evidence_multiplicity_to_log_t_scores = sort_t_score_dict(
			evidence_multiplicity_to_log_t_scores)

		fig = plt.figure(figsize=(9.5, 5))
		ax = fig.add_subplot(111)

		i = 0
		all_yticks = []
		labels = []

		for log_t_score_dict in [
			evidence_multiplicity_to_log_t_scores,
			evidence_code_to_log_t_scores,
			all_operons_log_t_scores,
			]:
			keys = [key for key in log_t_score_dict.keys()]
			log_t_scores = [value for value in log_t_score_dict.values()]
			n_samples = [len(t) for t in log_t_scores]
			yticks = np.arange(i, i + len(keys))
			all_yticks.extend(yticks)
			labels.extend([
				f'{EVIDENCE_CODE_TO_DESCRIPTIONS.get(key, key)} (n={n_sample})'
				for (key, n_sample) in zip(keys, n_samples)
				])

			for j, t_scores in enumerate(log_t_scores):
				parts = ax.violinplot(
					t_scores, positions=[i + j], widths=0.75,
					showextrema=False, vert=False)

				for pc in parts['bodies']:
					pc.set_facecolor('#cccccc')
					pc.set_edgecolor('none')
					pc.set_alpha(1)

				quartile1, median, quartile3 = np.percentile(
					t_scores, [25, 50, 75])

				ax.scatter(median, i + j, marker='o', color='white', s=10, zorder=3)
				ax.hlines(i + j, quartile1, quartile3, color='k', linestyle='-', lw=3)

			i = i + len(keys) + 1

		ax.set_xlim([0, 2])
		ax.set_ylim([-1, i - 1])
		ax.set_xticks([0, 0.5, 1, 1.5, 2])
		ax.xaxis.tick_top()
		ax.xaxis.set_label_position('top')
		ax.set_xlabel('log10(max t-score + 1)')
		ax.set_yticks(all_yticks)
		ax.set_yticklabels(labels)
		ax.spines["bottom"].set_visible(False)
		ax.spines["right"].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_evidence_codes', metadata)
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
