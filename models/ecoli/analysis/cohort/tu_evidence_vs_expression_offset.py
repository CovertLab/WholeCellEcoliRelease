"""
Generates a plot that compares how the distribution of p-values that quantify
how much the expression of each mRNA is offset from its expected values due to
the addition of transcription units versus the evidence code given for the
transcription units that include the given mRNA.
"""

from itertools import chain
import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import scipy.stats as st

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader


FIGSIZE = (9, 6)
NUMERICAL_ZERO = 1e-10
P_VALUE_MINIMUM = 1e-5
P_VALUE_CUTOFF = 1e-3
IGNORED_EVIDENCE_CODES = ['EV-ND']

EXPECTED_COUNT_CONDITION = 'basal'


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		simOutDir = os.path.join(cell_paths[0], 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Get mask for mRNA genes that are
		# i) Does not encode for ribosomal proteins or RNAPs
		# ii) not affected by manual overexpression/underexpression
		# to isolate the effects of operons on expression levels.
		is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		assert np.all(
			mRNA_cistron_ids == sim_data.process.transcription.cistron_data['id'][is_mRNA])

		mRNA_is_rnap_or_rprotein = np.logical_or(
				sim_data.process.transcription.cistron_data['is_RNAP'],
				sim_data.process.transcription.cistron_data['is_ribosomal_protein'])[is_mRNA]

		is_adjusted = np.zeros_like(is_mRNA, dtype=bool)
		all_rna_ids = sim_data.process.transcription.rna_data['id']
		for adjusted_cistron_id in sim_data.adjustments.rna_expression_adjustments.keys():
			# Include cistrons whose expression is adjusted because they belong
			# to the same TU as the cistron that is bumped up
			adjusted_rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(
				adjusted_cistron_id)
			adjusted_cistron_indexes = []
			for adjusted_rna_index in adjusted_rna_indexes:
				adjusted_cistron_indexes.extend(
					sim_data.process.transcription.rna_id_to_cistron_indexes(
						all_rna_ids[adjusted_rna_index]))

			is_adjusted[adjusted_cistron_indexes] = True

		mRNA_is_adjusted = is_adjusted[is_mRNA]
		mRNA_mask = np.logical_and(~mRNA_is_rnap_or_rprotein, ~mRNA_is_adjusted)

		# Get list of cistron IDs that satisfy these conditions
		plotted_mRNA_cistrons = []
		for i in np.where(mRNA_mask)[0]:
			plotted_mRNA_cistrons.append(mRNA_cistron_ids[i])

		# Get expected counts
		expected_counts = sim_data.process.transcription.cistron_expression[EXPECTED_COUNT_CONDITION][
			is_mRNA][mRNA_mask]

		# Normalize counts
		expected_counts /= expected_counts.sum()

		# Read actual counts
		all_actual_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts')

		# Get average count across all timesteps over all sims
		actual_counts = all_actual_counts[:, mRNA_mask].mean(axis=0)
		n_timesteps = all_actual_counts.shape[0]

		# Normalize counts
		actual_counts /= actual_counts.sum()

		# Get p-values for each cistron, assuming a Poissonian distribution
		z_scores = np.sqrt(n_timesteps/(expected_counts + NUMERICAL_ZERO)) * np.abs(
			actual_counts - expected_counts)
		p_values = st.norm.sf(z_scores)
		p_values[p_values < P_VALUE_MINIMUM] = P_VALUE_MINIMUM

		# Map each plotted cistron to list of evidence codes
		rna_ids = sim_data.process.transcription.rna_data['id']
		cistron_index_to_evidence_codes = {}
		for (i, cistron_id) in enumerate(plotted_mRNA_cistrons):
			if expected_counts[i] == 0 or actual_counts[i] == 0:
				continue
			evidence_codes = []
			rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id)
			for rna_index in rna_indexes:
				evidence_codes.extend(
					sim_data.process.transcription.rna_id_to_evidence_codes.get(rna_ids[rna_index][:-3], []))

			cistron_index_to_evidence_codes[i] = list(set(evidence_codes))

		# Map each evidence code and superclass to p-values of cistrons
		evidence_code_to_p_values = {}
		evidence_superclass_to_p_values = {}
		for (cistron_index, evidence_codes) in cistron_index_to_evidence_codes.items():
			for evidence_code in evidence_codes:
				if evidence_code in IGNORED_EVIDENCE_CODES:
					continue
				evidence_code_to_p_values.setdefault(
					evidence_code, []).append(p_values[cistron_index])
				evidence_superclass_to_p_values.setdefault(
					evidence_code.split('-')[1], []).append(p_values[cistron_index])

		# Sort dictionary keys by order of increasing cutoff rates
		evidence_code_to_p_values = dict(sorted(
			evidence_code_to_p_values.items(),
			key=lambda item: (np.array(item[1]) < P_VALUE_CUTOFF).sum()/len(item[1])
			))
		evidence_superclass_to_p_values = dict(sorted(
			evidence_superclass_to_p_values.items(),
			key=lambda item: (np.array(item[1]) < P_VALUE_CUTOFF).sum()/len(item[1])
			))

		plt.figure(figsize=FIGSIZE)

		x_labels = []
		for (i, (key, p_values)) in enumerate(chain(
				evidence_superclass_to_p_values.items(),
				evidence_code_to_p_values.items())):
			x_labels.append(key)
			x_jitter = np.random.normal(i, 0.05, len(p_values))
			cutoff_mask = np.array(p_values) < P_VALUE_CUTOFF
			plt.scatter(
				x_jitter[~cutoff_mask], np.array(p_values)[~cutoff_mask],
				s=1, c='k', alpha=0.5)
			plt.scatter(
				x_jitter[cutoff_mask], np.array(p_values)[cutoff_mask],
				s=1, c='b', alpha=0.5)

		plt.title('TU evidence codes vs p-values')
		plt.xticks(np.arange(len(x_labels)), x_labels, rotation=90)
		plt.ylabel('p-value')
		plt.yscale('log')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')




if __name__ == '__main__':
	Plot().cli()
