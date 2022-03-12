"""
Generates a scatter plot of mRNA cistron copy numbers that are expected from the
expression levels calculated in the ParCa vs actual copies in the simulations
for a set of simulations run with operons, and a set run without operons.
"""
import pickle
import os
from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import scipy.stats as st

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
P_VALUE_THRESHOLD = 1e-7

SEED = 0
EXPECTED_COUNT_CONDITION = 'basal'


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
		# iii) not monocistronic
		# to isolate the effects of operons on expression levels.
		is_mRNA = sim_data2.process.transcription.cistron_data['is_mRNA']
		assert np.all(
			mRNA_cistron_ids == sim_data2.process.transcription.cistron_data['id'][is_mRNA])

		mRNA_is_rnap_or_rprotein = np.logical_or(
				sim_data2.process.transcription.cistron_data['is_RNAP'],
				sim_data2.process.transcription.cistron_data['is_ribosomal_protein'])[is_mRNA]

		is_adjusted = np.zeros_like(is_mRNA, dtype=bool)
		all_rna_ids = sim_data2.process.transcription.rna_data['id']
		for adjusted_cistron_id in sim_data2.adjustments.rna_expression_adjustments.keys():
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

		polycistronic_cistron_indexes = []
		for rna_id in sim_data2.process.transcription.rna_data['id']:
			cistron_indexes = sim_data2.process.transcription.rna_id_to_cistron_indexes(rna_id)
			if len(cistron_indexes) > 1:
				polycistronic_cistron_indexes.extend(cistron_indexes)
		is_polycistronic = np.zeros(len(sim_data2.process.transcription.cistron_data), bool)
		if len(polycistronic_cistron_indexes) > 0:
			is_polycistronic[np.array(
				list(set(polycistronic_cistron_indexes)))] = True

		mRNA_is_polycistronic = is_polycistronic[is_mRNA]

		mask = np.logical_and(~mRNA_is_rnap_or_rprotein, ~mRNA_is_adjusted)
		mRNA_mask_polycistronic = mRNA_is_polycistronic[mask]

		# Get list of polycistron IDs that are plotted
		plotted_mRNA_ids = []
		for i in np.where(mask)[0]:
			plotted_mRNA_ids.append(mRNA_cistron_ids[i])

		fig = plt.figure(figsize=FIGSIZE)

		def read_data(ap):
			cell_paths = ap.get_cells(generation=[2, 3, 4, 5, 6, 7])
			all_actual_counts = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_cistron_counts')[:, mask]

			n_samples = 48
			n_bootstrap = 1000

			# Calculate mean from select timesteps
			counts = all_actual_counts[np.random.choice(all_actual_counts.shape[0], n_samples), :].mean(axis=0)

			# Estimate variance by bootstrapping
			sampled_counts = all_actual_counts[np.random.choice(all_actual_counts.shape[0], n_samples * n_bootstrap), :].reshape((n_bootstrap, n_samples, -1)).mean(axis=1)
			var = sampled_counts.var(axis=0)

			return counts, var

		counts1, v1 = read_data(ap1)
		counts2, v2 = read_data(ap2)

		# Normalize counts
		r = counts1.sum()/counts2.sum()
		counts2 *= r
		v2 *= r

		# Get statistical significance boundaries assuming a Poissonian
		# distribution
		z_score_threshold = st.norm.ppf(1 - P_VALUE_THRESHOLD/2)
		ub = counts1 + z_score_threshold*np.sqrt(v1 + v2)
		lb = counts1 - z_score_threshold*np.sqrt(v1 + v2)

		# Get mask for outliers
		outlier_mask = np.logical_or(counts2 > ub, counts2 < lb)

		ax = fig.add_subplot(1, 2, 1)
		ax.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
		ax.scatter(
			np.log10(counts1[np.logical_and(mRNA_mask_polycistronic, ~outlier_mask)] + 1),
			np.log10(counts2[np.logical_and(mRNA_mask_polycistronic, ~outlier_mask)] + 1),
			c='#cccccc', s=2,
			label=f'p ≥ {P_VALUE_THRESHOLD:g} (n = {np.logical_and(mRNA_mask_polycistronic, ~outlier_mask).sum():d})',
			clip_on=False)
		# Highlight outliers
		ax.scatter(
			np.log10(counts1[np.logical_and(mRNA_mask_polycistronic, outlier_mask)] + 1),
			np.log10(counts2[np.logical_and(mRNA_mask_polycistronic, outlier_mask)] + 1),
			c='b', s=2,
			label=f'p < {P_VALUE_THRESHOLD:g} (n = {np.logical_and(mRNA_mask_polycistronic, outlier_mask).sum():d})',
			clip_on=False)

		ax.set_title('Polycistronic genes')
		ax.set_xlabel('$\log_{10}$(mRNA copies + 1), old sims')
		ax.set_ylabel('$\log_{10}$(mRNA copies + 1), new sims')
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_position(("outward", 20))
		ax.spines["left"].set_position(("outward", 20))
		ax.set_xlim(BOUNDS)
		ax.set_ylim(BOUNDS)
		ax.legend(loc=2)

		ax = fig.add_subplot(1, 2, 2)
		ax.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
		ax.scatter(
			np.log10(counts1[np.logical_and(~mRNA_mask_polycistronic, ~outlier_mask)] + 1),
			np.log10(counts2[np.logical_and(~mRNA_mask_polycistronic, ~outlier_mask)] + 1),
			c='#cccccc', s=2,
			label=f'p ≥ {P_VALUE_THRESHOLD:g} (n = {np.logical_and(~mRNA_mask_polycistronic, ~outlier_mask).sum():d})',
			clip_on=False)
		# Highlight outliers
		ax.scatter(
			np.log10(counts1[np.logical_and(~mRNA_mask_polycistronic, outlier_mask)] + 1),
			np.log10(counts2[np.logical_and(~mRNA_mask_polycistronic, outlier_mask)] + 1),
			c='b', s=2,
			label=f'p < {P_VALUE_THRESHOLD:g} (n = {np.logical_and(~mRNA_mask_polycistronic, outlier_mask).sum():d})',
			clip_on=False)
		ax.set_title('Monocistronic genes')
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


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == '__main__':
	Plot().cli()
