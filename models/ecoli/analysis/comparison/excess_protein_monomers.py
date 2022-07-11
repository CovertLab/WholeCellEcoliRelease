"""
Comparison plot to compare the coefficients of variation of the expression
levels of protein monomers that constitute the same protein complex for sims
with/without polycistronic operons.
"""

from functools import reduce
import os
from typing import Tuple

import csv
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
# noinspection PyUnresolvedReferences
import numpy as np
from scipy import stats

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure,
    read_stacked_columns, read_stacked_bulk_molecules)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader



class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		# Load from sim_data
		all_subunit_ids = sim_data1.process.complexation.molecule_names
		all_complex_ids = sim_data1.process.complexation.ids_complexes
		stoich_matrix_monomers = sim_data1.process.complexation.stoich_matrix_monomers()
		monomer_id_to_cistron_id = {
			monomer['id']: monomer['cistron_id']
			for monomer in sim_data1.process.translation.monomer_data
			}
		operons = sim_data2.process.transcription.operons
		cistron_id_to_cistron_index = {
			cistron['id']: i
			for (i, cistron) in enumerate(sim_data2.process.transcription.cistron_data)
			}
		subunit_is_shared = ((stoich_matrix_monomers < 0).sum(axis=1) > 1)

		# Get ordering of monomers in MonomerCounts table
		monomer_counts_reader = TableReader(
			os.path.join(ap1.get_cells()[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)}

		def is_in_same_operon(subunit_ids):
			"""
			Returns True if the cistrons that each encode for the list of
			protein monomer subunits are in the same operon.
			"""
			# Get indexes of cistrons that encode for each monomer
			try:
				cistron_indexes = [
					cistron_id_to_cistron_index[monomer_id_to_cistron_id[monomer_id]]
					for monomer_id in subunit_ids]
			except KeyError:
				return False

			# Get list of operons that these cistrons belong to
			subunit_operons = [
				operon for operon in operons
				if np.any(np.isin(cistron_indexes, operon[0]))]

			# Return true if only one operon was found
			return len(subunit_operons) == 1

		complex_id_to_subunit_ids = {}
		complex_id_to_subunit_stoichs = {}
		complex_id_to_is_in_same_operon = {}
		complex_id_to_has_shared_subunits = {}

		# Get IDs and stoichiometries of subunits that compose each complex,
		# and whether those subunits are cotranscribed as part of a same operon
		for i, complex_id in enumerate(all_complex_ids):
			subunit_mask = stoich_matrix_monomers[:, i] < 0
			subunit_indexes = np.where(subunit_mask)[0]
			complex_id_to_subunit_stoichs[complex_id] = -stoich_matrix_monomers[:, i][subunit_mask]
			subunit_ids = [all_subunit_ids[i] for i in subunit_indexes]
			complex_id_to_subunit_ids[complex_id] = subunit_ids
			complex_id_to_is_in_same_operon[complex_id] = is_in_same_operon(
				subunit_ids)
			complex_id_to_has_shared_subunits[complex_id] = np.any(
				subunit_is_shared[subunit_indexes])

		# Get whether the stoichiometries of each complexes are unknown
		complex_id_to_stoich_unknown = {
			complex_id: is_unknown for (complex_id, is_unknown)
			in zip(all_complex_ids, sim_data1.process.complexation.reaction_stoichiometry_unknown)
		}

		def read_sims(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			all_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True)
			(all_complex_counts, ) = read_stacked_bulk_molecules(
				cell_paths, (all_complex_ids, ), ignore_exception=True)

			# Take averages across all sims and timepoints
			all_monomer_counts_mean = all_monomer_counts.mean(axis=0)
			all_complex_counts_mean = all_complex_counts.mean(axis=0)

			complex_id_to_cov = {}
			excess_monomer_counts = np.zeros_like(all_complex_counts_mean)

			# Loop through each protein complex
			for i, complex_id in enumerate(all_complex_ids):
				# Get stoichiometries and IDs of each monomer subunit
				subunit_stoichs = complex_id_to_subunit_stoichs[complex_id]
				subunit_ids = complex_id_to_subunit_ids[complex_id]

				# Skip homogeneous protein complexes
				if len(subunit_ids) == 1:
					continue

				# Get indexes of monomer subunits in the monomer counts table
				try:
					monomer_indexes = np.array([
						monomer_id_to_index[subunit_id] for subunit_id in subunit_ids
						])
				except KeyError:
					# Skip complexes with non-protein monomers or monomers that
					# were excluded from the model
					continue

				# Get counts of monomers and rescale with stoichiometries
				monomer_counts = all_monomer_counts_mean[monomer_indexes]
				rescaled_monomer_counts = monomer_counts / subunit_stoichs

				# Skip complexes with zero counts
				if all_complex_counts_mean[i] == 0:
					continue

				# Calculate coefficient of variation
				complex_id_to_cov[complex_id] = (
					rescaled_monomer_counts.std()
					/ rescaled_monomer_counts.mean())

				# Calculate excess monomer counts for this complex
				excess_monomer_counts[i] = monomer_counts.sum() - all_complex_counts_mean[i]*subunit_stoichs.sum()

			return complex_id_to_cov, all_monomer_counts_mean, all_complex_counts_mean, excess_monomer_counts

		cov1, m1, c1, em1 = read_sims(ap1)
		cov2, m2, c2, em2 = read_sims(ap2)

		# Select complexes that have coefficients of variation calculated for
		# both sims and has cotranscribed subunits, known stoichiometries, and
		# no shared subunits with other complexes
		complexes_to_plot = [
			complex_id for complex_id in cov1.keys()
			if complex_id in cov2
				and complex_id_to_is_in_same_operon[complex_id]
				and not complex_id_to_stoich_unknown[complex_id]
				and not complex_id_to_has_shared_subunits[complex_id]
			]

		# Calculate total counts of relevant complexes and monomers
		complex_id_to_index = {
			complex_id: i for (i, complex_id) in enumerate(all_complex_ids)}
		complex_indexes = np.array([
			complex_id_to_index[complex_id] for complex_id in complexes_to_plot])
		monomer_indexes = []
		for complex_id in complexes_to_plot:
			monomer_indexes.extend([
				monomer_id_to_index[subunit_id] for subunit_id
				in complex_id_to_subunit_ids[complex_id]])
		monomer_indexes = np.array(monomer_indexes)

		m_total1 = m1[monomer_indexes].sum()
		m_total2 = m2[monomer_indexes].sum()
		c_total1 = c1[complex_indexes].sum()
		c_total2 = c2[complex_indexes].sum()
		em_total1 = em1[complex_indexes].sum()
		em_total2 = em2[complex_indexes].sum()

		# Plot comparison of coefficient distributions
		fig = plt.figure(figsize=(9.05, 4))
		gs = gridspec.GridSpec(
			2, 4, width_ratios=(4, 1, 4, 1.9), height_ratios=(1, 4))

		def draw_plot(p1, p2, grid_i, grid_j, y_max):
			scatter_ax = fig.add_subplot(gs[grid_i, grid_j])
			scatter_ax.plot([0, 2], [0, 2], ls='--', lw=1, c='k', alpha=0.1)
			scatter_ax.scatter(
				p1, p2,
				alpha=0.4, s=10, c='#555555', clip_on=False, edgecolors='none')

			scatter_ax.set_xlim([0, 2.0])
			scatter_ax.set_ylim([0, 2.0])
			scatter_ax.set_xticks([0, 2.0])
			scatter_ax.set_yticks([0, 2.0])
			scatter_ax.set_xlabel('Reference')
			scatter_ax.set_ylabel('Input')

			scatter_ax.spines["top"].set_visible(False)
			scatter_ax.spines["right"].set_visible(False)
			scatter_ax.spines["bottom"].set_position(("outward", 20))
			scatter_ax.spines["left"].set_position(("outward", 20))

			x = np.linspace(0, 2.0, 1000)
			kde1 = stats.gaussian_kde(p1)
			kde2 = stats.gaussian_kde(p2)

			hist1_ax = fig.add_subplot(gs[grid_i - 1, grid_j], sharex=scatter_ax)
			hist1_ax.fill_between(x, kde1(x), alpha=0.5)
			hist1_ax.axvline(p1.mean(), lw=2, ls='--', c='#555555')
			hist1_ax.set_xlim([0, 2.0])
			hist1_ax.set_ylim([0, y_max])
			hist1_ax.set_yticks([])
			hist1_ax.spines["top"].set_visible(False)
			hist1_ax.spines["right"].set_visible(False)
			hist1_ax.spines["left"].set_visible(False)
			hist1_ax.spines["bottom"].set_visible(False)
			plt.setp(hist1_ax.get_xaxis(), visible=False)

			hist2_ax = fig.add_subplot(gs[grid_i, grid_j + 1], sharey=scatter_ax)
			hist2_ax.fill_betweenx(x, kde2(x), fc='C1', alpha=0.5)
			hist2_ax.axhline(p2.mean(), lw=2, ls='--', c='#555555')
			hist2_ax.set_ylim([0, 2.0])
			hist2_ax.set_xlim([0, y_max])
			hist2_ax.set_xticks([])
			hist2_ax.spines["top"].set_visible(False)
			hist2_ax.spines["right"].set_visible(False)
			hist2_ax.spines["left"].set_visible(False)
			hist2_ax.spines["bottom"].set_visible(False)
			plt.setp(hist2_ax.get_yaxis(), visible=False)

		y1 = np.array([
			cov1[complex_id] for complex_id in complexes_to_plot])
		y2 = np.array([
			cov2[complex_id] for complex_id in complexes_to_plot])
		draw_plot(y1, y2, 1, 0, 3)

		# Plot comparisons of total monomer and complex counts
		ax_monomer = fig.add_subplot(gs[0:2, 2])
		ax_monomer.bar(-0.15, m_total1, width=0.3, alpha=0.5, label='reference')
		ax_monomer.bar(0.15, m_total2, width=0.3, alpha=0.5, label='input')
		ax_monomer.bar(0.85, em_total1, width=0.3, alpha=0.5, color='C0')
		ax_monomer.bar(1.15, em_total2, width=0.3, alpha=0.5, color='C1')
		ax_monomer.set_xticks([0, 1])
		ax_monomer.set_xticklabels(['all\nsubunits', 'excess\nsubunits'])
		ax_monomer.set_xlim([-0.5, 1.5])
		ax_monomer.set_ylabel('Average total counts')
		ax_monomer.set_ylim([0, 250000])
		ax_monomer.spines['top'].set_visible(False)
		ax_monomer.spines['right'].set_visible(False)
		plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

		ax_complex = fig.add_subplot(gs[0:2, 3])
		ax_complex.bar(-0.15, c_total1, width=0.3, alpha=0.5, label='reference')
		ax_complex.bar(0.15, c_total2, width=0.3, alpha=0.5, label='input')
		ax_complex.set_xticks([0])
		ax_complex.set_xticklabels(['complexes'])
		ax_complex.set_xlim([-0.5, 0.5])
		ax_complex.set_ylabel('Average counts')
		ax_complex.set_ylim([0, 20000])
		ax_complex.set_yticks([0, 10000, 20000])
		ax_complex.spines['top'].set_visible(False)
		ax_complex.spines['right'].set_visible(False)
		plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

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
