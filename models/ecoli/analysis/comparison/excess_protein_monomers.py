"""
Comparison plot to compare the Gini coefficients of the expression levels of
protein monomers that constitute the same protein complex for sims with/without
polycistronic operons.
"""

from functools import reduce
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
from wholecell.analysis.analysis_tools import (exportFigure,
    read_stacked_columns, read_stacked_bulk_molecules)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


SEED = 0

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		np.random.seed(SEED)

		# Load from sim_data
		all_subunit_ids = sim_data1.process.complexation.molecule_names
		all_complex_ids = sim_data1.process.complexation.ids_complexes
		stoich_matrix_monomers = sim_data1.process.complexation.stoich_matrix_monomers()
		monomer_id_to_cistron_id = {
			monomer['id']: monomer['cistron_id']
			for monomer in sim_data1.process.translation.monomer_data
			}
		cistron_id_to_rna_indexes = sim_data2.process.transcription.cistron_id_to_rna_indexes

		# Get ordering of monomers in MonomerCounts table
		monomer_counts_reader = TableReader(
			os.path.join(ap1.get_cells()[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)}

		def is_cotranscribed(subunit_ids):
			"""
			Returns True if the cistrons that each encode for the list of
			protein monomer subunits are ever found on the same transcription
			unit.
			"""
			# Get IDs of cistrons that encode for each monomer
			try:
				cistron_ids = [
					monomer_id_to_cistron_id[monomer_id]
					for monomer_id in subunit_ids]
			except KeyError:
				return False

			# Get indexes of transcription units that contain each cistron
			rna_indexes = [
				cistron_id_to_rna_indexes(cistron_id)
				for cistron_id in cistron_ids]

			# Find any overlapping transcription unit indexes
			return len(reduce(np.intersect1d, rna_indexes)) > 0

		complex_id_to_subunit_ids = {}
		complex_id_to_subunit_stoichs = {}
		complex_id_to_is_cotranscribed = {}

		# Get IDs and stoichiometries of subunits that compose each complex,
		# and whether those subunits are cotranscribed as part of a same operon
		for i, complex_id in enumerate(all_complex_ids):
			subunit_mask = stoich_matrix_monomers[:, i] < 0
			subunit_indexes = np.where(subunit_mask)[0]
			complex_id_to_subunit_stoichs[complex_id] = -stoich_matrix_monomers[:, i][subunit_mask]
			subunit_ids = [all_subunit_ids[i] for i in subunit_indexes]
			complex_id_to_subunit_ids[complex_id] = subunit_ids
			complex_id_to_is_cotranscribed[complex_id] = is_cotranscribed(
				subunit_ids)

		# Get whether the stoichiometries of each complexes are unknown
		complex_id_to_stoich_unknown = {
			complex_id: is_unknown for (complex_id, is_unknown)
			in zip(all_complex_ids, sim_data1.process.complexation.reaction_stoichiometry_unknown)
		}

		def read_sims(ap):
			all_monomer_counts = read_stacked_columns(
				ap.get_cells(), 'MonomerCounts', 'monomerCounts')
			(all_complex_counts, ) = read_stacked_bulk_molecules(
				ap.get_cells(), (all_complex_ids, ))

			# Take averages across all sims and timepoints
			all_monomer_counts_mean = all_monomer_counts.mean(axis=0)
			all_complex_counts_mean = all_complex_counts.mean(axis=0)

			complex_id_to_gini_coeff = {}

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
				if np.all(rescaled_monomer_counts == 0):
					continue

				# Calculate Gini coefficient
				complex_id_to_gini_coeff[complex_id] = (
					np.sum(np.abs(np.expand_dims(rescaled_monomer_counts, 0)
					    - np.expand_dims(rescaled_monomer_counts, 1)))
					/ (2 * len(rescaled_monomer_counts) * rescaled_monomer_counts.sum())
				)

			return complex_id_to_gini_coeff, all_monomer_counts_mean, all_complex_counts_mean

		gini_coeff1, m1, c1 = read_sims(ap1)
		gini_coeff2, m2, c2 = read_sims(ap2)

		# Select complexes that have Gini coefficients calculated for both sims
		# and has cotranscribed subunits and known stoichiometries
		complexes_to_plot = [
			complex_id for complex_id in gini_coeff1.keys()
			if complex_id in gini_coeff2
			   and complex_id_to_is_cotranscribed[complex_id]
			   and not complex_id_to_stoich_unknown[complex_id]
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

		# Plot comparison of Gini coefficient distributions
		plt.figure(figsize=(6, 4))
		gs = gridspec.GridSpec(1, 2)

		ax_gini = plt.subplot(gs[0, 0])
		y1 = np.array([
			gini_coeff1[complex_id] for complex_id in complexes_to_plot])
		y2 = np.array([
			gini_coeff2[complex_id] for complex_id in complexes_to_plot])

		x_jitter1 = np.random.normal(0, 0.02, len(complexes_to_plot))
		x_jitter2 = np.random.normal(0.8, 0.02, len(complexes_to_plot))

		ax_gini.scatter(x_jitter1, y1, s=5)
		ax_gini.scatter(x_jitter2, y2, s=5)

		ax_gini.violinplot(y1, positions=[0], showextrema=False)
		ax_gini.violinplot(y2, positions=[0.8], showextrema=False)

		ax_gini.set_xticks([0, 0.8])
		ax_gini.set_xticklabels(['reference', 'input'])
		ax_gini.set_yticks([0, 1])
		ax_gini.set_xlim([-0.5, 1.3])
		ax_gini.set_ylim([0, 1.1])
		ax_gini.set_ylabel('Gini coefficient')

		ax_gini.spines['top'].set_visible(False)
		ax_gini.spines['right'].set_visible(False)

		# Plot comparisons of total monomer and complex counts
		ax_monomer = plt.subplot(gs[0, 1])
		ax_monomer.bar(-0.15, m_total1, width=0.3)
		ax_monomer.bar(0.15, m_total2, width=0.3)
		ax_monomer.set_xticks([0, 1])
		ax_monomer.set_xticklabels(['subunits', 'complexes'])
		ax_monomer.set_xlim([-0.5, 1.5])
		ax_monomer.set_ylabel('Total subunit counts')
		ax_monomer.set_ylim([0, 600000])
		ax_monomer.spines['top'].set_visible(False)
		plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

		ax_complex = ax_monomer.twinx()
		ax_complex.bar(0.85, c_total1, width=0.3, label='reference')
		ax_complex.bar(1.15, c_total2, width=0.3, label='input')
		ax_complex.set_ylabel('Total complex counts')
		ax_complex.set_ylim([0, 50000])
		ax_complex.spines['top'].set_visible(False)
		plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
		ax_complex.legend(loc=1, prop={'size': 8})

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
