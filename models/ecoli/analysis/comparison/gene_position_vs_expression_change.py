"""
Plots the relationship between the position of a gene within an operon and the
changes in its expression levels between simulations with and without operons.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os
import re
from scipy.stats import pearsonr

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts
from wholecell.utils import units

NUMERICAL_ZERO = 1e-10

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.

		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		transcription = sim_data2.process.transcription
		operons = transcription.operons
		cistron_start_pos = transcription.cistron_data['replication_coordinate']
		cistron_is_forward = transcription.cistron_data['is_forward']
		cistron_lengths = transcription.cistron_data['length'].asNumber(units.nt)
		cistron_end_pos = cistron_start_pos + 2*(cistron_is_forward - 0.5)*cistron_lengths

		# Get genomic coordinates of the center points of each operon
		operon_positions = []
		polycistronic_indexes = []
		for operon in operons:
			cistron_indexes = operon[0]
			boundaries = np.concatenate((
				cistron_start_pos[cistron_indexes],
				cistron_end_pos[cistron_indexes]
				))
			operon_positions.append((boundaries.max() + boundaries.min())/2)

			if len(cistron_indexes) > 1:
				polycistronic_indexes.extend(cistron_indexes)

		cistron_index_to_operon_pos = {}
		for operon, pos in zip(operons, operon_positions):
			for cistron_index in operon[0]:
				cistron_index_to_operon_pos[cistron_index] = pos

		# Calculate offsets between gene start site and operon centers
		cistron_id_to_index = {
			cistron_id: i for i, cistron_id in enumerate(transcription.cistron_data['id'])
			}
		cistron_indexes_reordered = np.array([
			cistron_id_to_index[cistron_id]
			for cistron_id in sim_data2.process.translation.monomer_data['cistron_id']
			])
		monomer_is_polycistronic = np.array([
			(cistron_index in polycistronic_indexes)
			for cistron_index in cistron_indexes_reordered
			])
		monomer_start_pos = cistron_start_pos[cistron_indexes_reordered]
		monomer_is_forward = cistron_is_forward[cistron_indexes_reordered]
		monomer_operon_pos = np.array([
			cistron_index_to_operon_pos[i] for i in cistron_indexes_reordered
			])
		pos_diff = 2 * (monomer_is_forward - 0.5) * (
				monomer_start_pos - monomer_operon_pos)

		def read_sim_protein_counts(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True).mean(axis=0)

			return monomer_counts

		# Read monomer counts
		monomer_counts1 = read_sim_protein_counts(ap1)
		monomer_counts2 = read_sim_protein_counts(ap2)

		count_ratios = (monomer_counts2 + NUMERICAL_ZERO) / (monomer_counts1 + NUMERICAL_ZERO)
		is_expressed = np.logical_and(
			monomer_counts1 > 10, monomer_counts2 > 10)

		plt.figure(figsize=(4, 4.2))
		ax = plt.subplot(1, 1, 1)

		ax.scatter(
			pos_diff[np.logical_and(is_expressed, monomer_is_polycistronic)],
			count_ratios[np.logical_and(is_expressed, monomer_is_polycistronic)],
			alpha=0.4, s=10, c='#555555', clip_on=False, edgecolors='none')
		ax.set_xlim([-10000, 10000])
		ax.set_ylim([1e-2, 1e2])
		ax.set_yscale('log')
		ax.set_xlabel('Position from center of operon (nt)')
		ax.set_ylabel('Fold change in protein counts with operons')
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_position(("outward", 20))
		ax.spines["left"].set_position(("outward", 20))

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
