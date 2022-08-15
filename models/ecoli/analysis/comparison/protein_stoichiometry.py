"""
Comparison plot to compare how far the stoichiometry of proteins that are
encoded from the same operon deviates from their mean values for sims
with/without polycistronic operons.
"""

import csv
import os
from typing import Tuple

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
# noinspection PyUnresolvedReferences
import numpy as np
from scipy import stats

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader

TU_ID = 'TU0-1123[c]'


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		# Load from sim_data
		operon_cistron_ids = sim_data2.process.transcription.cistron_data['id']
		cistron_id_to_monomer_id = {
			monomer['cistron_id']: monomer['id'] for monomer
			in sim_data2.process.translation.monomer_data
			}
		operons = sim_data2.process.transcription.operons

		# Load monomer ID attribute
		monomer_counts_reader = TableReader(
			os.path.join(ap1.get_cells()[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)}

		# Load list of tuples of monomer indexes that are in the same
		# polycistronic operon
		operon_monomer_indexes = []

		for operon in operons:
			cistron_indexes = operon[0]
			if len(cistron_indexes) < 2:
				continue
			operon_monomer_indexes.append(
				np.array([monomer_id_to_index[cistron_id_to_monomer_id[operon_cistron_ids[i]]] for i in cistron_indexes])
				)

		def read_cv(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			all_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True)

			# Take averages across all sims and timepoints
			all_monomer_counts_mean = all_monomer_counts.mean(axis=0)

			# Initialize array of maximum coefficients of variation
			is_constitutive = np.zeros(len(operon_monomer_indexes), dtype=np.bool)
			max_cv = np.zeros(len(operon_monomer_indexes))

			# Loop through each operon
			for i, monomer_indexes in enumerate(operon_monomer_indexes):
				# Get counts of each monomer subunit
				monomer_counts = all_monomer_counts[:, monomer_indexes]
				avg_monomer_counts = all_monomer_counts_mean[monomer_indexes]

				if np.any(monomer_counts == 0):
					continue

				# Rescale counts with average counts
				rescaled_monomer_counts = monomer_counts / avg_monomer_counts

				# Calculate coefficients of variation at every timestep
				cv = rescaled_monomer_counts.std(axis=1)/rescaled_monomer_counts.mean(axis=1)

				is_constitutive[i] = True
				max_cv[i] = cv.max()

			return max_cv, is_constitutive

		max_cv1, is_constitutive1 = read_cv(ap1)
		max_cv2, is_constitutive2 = read_cv(ap2)

		# Plot comparison of coefficient distributions
		fig = plt.figure(figsize=(4.5, 4.4))
		gs = gridspec.GridSpec(
			2, 2, width_ratios=(4, 1), height_ratios=(1, 4))

		def draw_scatter_plot(p1, p2, grid_i, grid_j, y_max):
			scatter_ax = fig.add_subplot(gs[grid_i, grid_j])
			scatter_ax.plot([0, 1.5], [0, 1.5], ls='--', lw=1, c='k', alpha=0.1)
			scatter_ax.scatter(
				p1, p2,
				alpha=0.4, s=10, c='#555555', clip_on=False, edgecolors='none')

			scatter_ax.set_xlim([0, 1.5])
			scatter_ax.set_ylim([0, 1.5])
			scatter_ax.set_xticks([0, 1.5])
			scatter_ax.set_yticks([0, 1.5])
			scatter_ax.set_xlabel('Reference')
			scatter_ax.set_ylabel('Input')

			scatter_ax.spines["top"].set_visible(False)
			scatter_ax.spines["right"].set_visible(False)
			scatter_ax.spines["bottom"].set_position(("outward", 20))
			scatter_ax.spines["left"].set_position(("outward", 20))

			x = np.linspace(0, 1.5, 1000)
			kde1 = stats.gaussian_kde(p1)
			kde2 = stats.gaussian_kde(p2)

			hist1_ax = fig.add_subplot(gs[grid_i - 1, grid_j], sharex=scatter_ax)
			hist1_ax.fill_between(x, kde1(x), alpha=0.5)
			hist1_ax.axvline(p1.mean(), lw=2, ls='--', c='#555555')
			hist1_ax.set_xlim([0, 1.5])
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
			hist2_ax.set_ylim([0, 1.5])
			hist2_ax.set_xlim([0, y_max])
			hist2_ax.set_xticks([])
			hist2_ax.spines["top"].set_visible(False)
			hist2_ax.spines["right"].set_visible(False)
			hist2_ax.spines["left"].set_visible(False)
			hist2_ax.spines["bottom"].set_visible(False)
			plt.setp(hist2_ax.get_yaxis(), visible=False)

		plotted_operons = np.logical_and(is_constitutive1, is_constitutive2)
		draw_scatter_plot(
			max_cv1[plotted_operons], max_cv2[plotted_operons],
			1, 0, 3)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		operon_cistron_ids = [
			sim_data2.process.transcription.cistron_data['id'][i]
			for i in sim_data2.process.transcription.rna_id_to_cistron_indexes(TU_ID)]
		cistron_monomer_ids = [
			cistron_id_to_monomer_id[cistron_id] for cistron_id in operon_cistron_ids]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id']
			for cistron in sim_data2.process.transcription.cistron_data
			}
		gene_id_to_gene_symbol = {
			gene['name']: gene['symbol']
			for gene in sim_data2.process.replication.gene_data
			}
		gene_symbols = [
			gene_id_to_gene_symbol[cistron_id_to_gene_id[cistron_id]]
			for cistron_id in operon_cistron_ids]

		def read_monomer_trace(ap):
			cell_paths = ap.get_cells(seed=[14])

			simOutDir = os.path.join(cell_paths[0], "simOut")
			monomer_counts_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))
			all_monomer_ids = monomer_counts_reader.readAttribute('monomerIds')

			# Get indexes of monomers translated from the TU
			monomer_indexes = [
				all_monomer_ids.index(monomer_id) for monomer_id in cistron_monomer_ids]

			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time')
			gen_start_time = read_stacked_columns(
				cell_paths, 'Main', 'time', fun=lambda x: x[0])
			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts')[:, monomer_indexes]
			monomer_counts_ratios = (monomer_counts[:, 1:].T / monomer_counts[:, 0]).T

			return time, gen_start_time, monomer_counts, monomer_counts_ratios

		t1, gen_start_time1, mc1, mcr1 = read_monomer_trace(ap1)
		t2, gen_start_time2, mc2, mcr2 = read_monomer_trace(ap2)

		plt.figure(figsize=(8, 6))

		def draw_trace_plot(c, cr, t, gst, ax_index):
			# Plot monomer counts for reference sims
			ax1 = plt.subplot(4, 1, ax_index)
			for i, gene_symbol in enumerate(gene_symbols):
				ax1.plot(t / 60, c[:, i], label=gene_symbol, clip_on=False)
			if ax_index == 1:
				ax1.legend(loc='center left', prop={'size': 8}, bbox_to_anchor=(1, 0.5))
			ax1.set_ylabel('Protein\ncounts')
			ax1.spines["top"].set_visible(False)
			ax1.spines["right"].set_visible(False)
			ax1.spines["bottom"].set_position(("outward", 10))
			ax1.spines["left"].set_position(("outward", 10))
			ax1.spines["bottom"].set_visible(False)
			ax1.get_xaxis().set_visible(False)
			ax1.set_xlim([0, t[-1] / 60])
			ax1.set_ylim([0, 800])
			ax1.set_yticks([0, 800])

			# Plot monomer count ratios for reference sims
			ax2 = plt.subplot(4, 1, ax_index + 1)
			for i, color in enumerate(['C1', 'C2', 'C3', 'C4']):
				ax2.plot(t / 60, cr[:, i], label=f'{gene_symbols[i + 1]}:{gene_symbols[0]}', clip_on=False, color=color, ls='--')
			if ax_index == 1:
				ax2.legend(loc='center left', prop={'size': 8}, bbox_to_anchor=(1, 0.5))
			ax2.set_ylabel('Ratios')
			ax2.spines["top"].set_visible(False)
			ax2.spines["right"].set_visible(False)
			ax2.spines["bottom"].set_position(("outward", 10))
			ax2.spines["left"].set_position(("outward", 10))
			ax2.set_xlim([0, t[-1] / 60])
			ax2.set_yscale('log')
			ax2.set_ylim([1, 100])
			ax2.set_yticks([1, 10, 100])
			ax2.set_xticks(list(gst / 60) + [t[-1] / 60])
			ax2.set_xticklabels(np.arange(len(gst) + 1))
			ax2.set_xlabel('Generations')

		draw_trace_plot(mc1, mcr1, t1, gen_start_time1, 1)
		draw_trace_plot(mc2, mcr2, t2, gen_start_time2, 3)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_count_trace', metadata)
		plt.close('all')

		# Write table
		gene_group_indexes = [
			operon[0] for operon in sim_data2.process.transcription.operons
			if len(operon[0]) > 1]
		gene_ids = sim_data1.process.transcription.cistron_data['gene_id']
		gene_id_to_gene_symbol = {
			gene['name']: gene['symbol']
			for gene in sim_data1.process.replication.gene_data
			}
		constituent_gene_names = [
			[gene_id_to_gene_symbol[gene_ids[i]] for i in gene_indexes]
			for gene_indexes in gene_group_indexes
			]

		with open(os.path.join(plotOutDir, plotOutFileName + '_table.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'constituent_gene_names', 'c.v. without operons',
				'c.v. with operons'
				])

			for plotted, gene_names, c1, c2 in zip(
					plotted_operons, constituent_gene_names, max_cv1, max_cv2
					):
				if plotted:
					writer.writerow([
						', '.join(gene_names), c1, c2
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
