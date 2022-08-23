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
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure,
    read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


DIFF_THRESHOLD = 0.05
WRITE_TABLE = False

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, inputDir1, plotOutDir, plotOutFileName, inputDir2, unused, metadata):
		# operons="off"
		ap1, sim_data1, _ = self.setup(inputDir1)
		# operons="on"
		ap2, sim_data2, _ = self.setup(inputDir2)

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

			# Take averages across all sims and timepoints
			all_monomer_counts_mean = all_monomer_counts.mean(axis=0)

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

			return complex_id_to_gini_coeff

		gini_coeff1 = read_sims(ap1)
		gini_coeff2 = read_sims(ap2)

		# Select complexes that have Gini coefficients calculated for both sims
		# and has cotranscribed subunits and known stoichiometries
		complexes_to_plot = [
			complex_id for complex_id in gini_coeff1.keys()
			if complex_id in gini_coeff2
			   and complex_id_to_is_cotranscribed[complex_id]
			   and not complex_id_to_stoich_unknown[complex_id]
			]
		complex_ids_increased_gini = []

		# Plot comparison
		plt.figure(figsize=(4, 6))

		for complex_id in complexes_to_plot:
			# Highlight complexes whose coefficients increased with operons
			if gini_coeff2[complex_id] - gini_coeff1[complex_id] > DIFF_THRESHOLD:
				c = '#222222'
				lw = 2
				zorder = 5
				complex_ids_increased_gini.append(complex_id)
			else:
				c = '#cccccc'
				lw = 1
				zorder = 0
			plt.plot(
				[0, 1], [gini_coeff1[complex_id], gini_coeff2[complex_id]],
				color=c, lw=lw, zorder=zorder)

		plt.xticks([0, 1], ["off", "on"])
		plt.yticks([0, 1])
		plt.xlim([-0.2, 1.2])
		plt.ylim([-0.1, 1.1])
		plt.ylabel('Gini coefficient')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# If True, generate table of relevant parameters for complexes whose
		# coefficients increased significantly with operons
		if WRITE_TABLE:
			cistron_exp_off = sim_data2.process.transcription.cistron_expression['basal']
			cistron_exp_on = sim_data2.process.transcription.cistron_tu_mapping_matrix.dot(
				sim_data2.process.transcription.rna_expression['basal'])
			cistron_id_to_index = {
				cistron_id: i for i, cistron_id
				in enumerate(sim_data2.process.transcription.cistron_data['id'])}
			translation_efficiencies = sim_data1.process.translation.translation_efficiencies_by_monomer

			with open(os.path.join(plotOutDir, plotOutFileName) + '.tsv', 'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow([
					'complex_id', 'gini_operons_off', 'gini_operons_on',
					'subunit_ids', 'subunit_stoichs', 'subunit_exp_off',
					'subunit_exp_on', 'subunit_translation_eff'
					])

				# Add one row for each complex
				for complex_id in complex_ids_increased_gini:
					gini_operons_off = gini_coeff1[complex_id]
					gini_operons_on = gini_coeff2[complex_id]
					subunit_ids = complex_id_to_subunit_ids[complex_id]
					subunit_stoichs = complex_id_to_subunit_stoichs[complex_id]

					subunit_cistron_ids = [
						monomer_id_to_cistron_id[monomer_id] for monomer_id in subunit_ids]
					cistron_indexes = np.array([
						cistron_id_to_index[cistron_id] for cistron_id in subunit_cistron_ids])
					exp_off = cistron_exp_off[cistron_indexes]
					exp_on = cistron_exp_on[cistron_indexes]

					subunit_indexes = np.array([
						monomer_id_to_index[monomer_id] for monomer_id in subunit_ids])
					translation_effs = translation_efficiencies[subunit_indexes]

					writer.writerow([
						f'{complex_id}', f'{gini_operons_off:.2f}',
						f'{gini_operons_on:.2f}',
						', '.join([subunit_id for subunit_id in subunit_ids]),
						', '.join([f'{int(stoich):d}' for stoich in subunit_stoichs]),
						', '.join([f'{e:.2g}' for e in exp_off]),
						', '.join([f'{e:.2g}' for e in exp_on]),
						', '.join([f'{e:.2g}' for e in translation_effs])
						])


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
