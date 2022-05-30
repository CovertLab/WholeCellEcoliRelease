"""
Template for comparison analysis plots
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader

TU_ID = 'TU0-45225[c]'


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off", input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.
		ap1, _, _ = self.setup(reference_sim_dir)
		ap2, sim_data, _ = self.setup(input_sim_dir)

		cistron_ids = [
			sim_data.process.transcription.cistron_data['id'][i]
			for i in sim_data.process.transcription.rna_id_to_cistron_indexes(TU_ID)]
		cistron_id_to_monomer_id = {
			monomer['cistron_id']: monomer['id']
			for monomer in sim_data.process.translation.monomer_data
			}
		monomer_ids = [
			cistron_id_to_monomer_id[cistron_id] for cistron_id in cistron_ids]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id']
			for cistron in sim_data.process.transcription.cistron_data
			}
		gene_id_to_gene_symbol = {
			gene['name']: gene['symbol']
			for gene in sim_data.process.replication.gene_data
			}
		gene_symbols = [
			gene_id_to_gene_symbol[cistron_id_to_gene_id[cistron_id]]
			for cistron_id in cistron_ids]

		def read_sims(ap):
			cell_paths = ap.get_cells(seed=[0])

			simOutDir = os.path.join(cell_paths[0], "simOut")
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
			mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')
			monomer_counts_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))
			all_monomer_ids = monomer_counts_reader.readAttribute('monomerIds')

			# Get indexes of cistrons that constitute the TU
			cistron_indexes = [
				mRNA_cistron_ids.index(cistron_id) for cistron_id in cistron_ids]

			# Get indexes of monomers translated from the TU
			monomer_indexes = [
				all_monomer_ids.index(monomer_id) for monomer_id in monomer_ids]

			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time')
			gen_start_time = read_stacked_columns(
				cell_paths, 'Main', 'time', fun=lambda x: x[0])
			cistron_counts = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_cistron_counts')[:, cistron_indexes]
			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts')[:, monomer_indexes]

			return time, gen_start_time, cistron_counts, monomer_counts

		time_off, gen_start_time_off, cistron_counts_off, monomer_counts_off = read_sims(ap1)
		time_on, gen_start_time_on, cistron_counts_on, monomer_counts_on = read_sims(ap2)

		plt.figure(figsize=(7.5, 3))

		# Plot counts of each cistron when operon="off"
		ax1 = plt.subplot(2, 2, 1)
		ax1.plot(time_off/60, cistron_counts_off, clip_on=False)
		ax1.set_ylabel('RNA Counts')
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 15))
		ax1.spines["left"].set_position(("outward", 15))
		ax1.spines["bottom"].set_visible(False)
		ax1.get_xaxis().set_visible(False)
		ax1.set_xlim([0, time_off[-1]/60])
		ax1.set_ylim([0, 2])
		ax1.set_yticks([0, 1, 2])

		# Plot counts of each cistron when operon="on"
		ax2 = plt.subplot(2, 2, 2, sharey=ax1)
		ax2.plot(
			time_on/60, cistron_counts_on[:, 0],
			color='#555555', label=TU_ID[:-3], clip_on=False)
		ax2.legend(loc=1, prop={'size': 8})
		ax2.spines["top"].set_visible(False)
		ax2.spines["right"].set_visible(False)
		ax2.spines["bottom"].set_position(("outward", 15))
		ax2.spines["left"].set_position(("outward", 15))
		ax2.spines["bottom"].set_visible(False)
		ax2.get_xaxis().set_visible(False)
		ax2.spines["left"].set_visible(False)
		ax2.get_yaxis().set_visible(False)
		ax2.set_xlim([0, time_on[-1] / 60])

		# Plot counts of each protein when operon="off"
		ax3 = plt.subplot(2, 2, 3)
		ax3.plot(time_off/60, monomer_counts_off, clip_on=False)
		ax3.set_ylabel('Protein Counts')
		ax3.set_xlabel('Generations')
		ax3.spines["top"].set_visible(False)
		ax3.spines["right"].set_visible(False)
		ax3.spines["bottom"].set_position(("outward", 15))
		ax3.spines["left"].set_position(("outward", 15))
		ax3.set_xlim([0, time_off[-1] / 60])
		ax3.set_ylim([0, 500])
		ax3.set_yticks([0, 500])
		ax3.set_xticks(list(gen_start_time_off / 60) + [time_off[-1]/60])
		ax3.set_xticklabels(np.arange(len(gen_start_time_off) + 1))

		# Plot counts of each protein when operon="on"
		ax4 = plt.subplot(2, 2, 4, sharey=ax3)
		for i, gene_id in enumerate(gene_symbols):
			ax4.plot(
				time_on/60, monomer_counts_on[:, i],
				label=gene_id, clip_on=False)
		ax4.legend(loc=1, prop={'size': 8})
		ax4.set_xlabel('Generations')
		ax4.spines["top"].set_visible(False)
		ax4.spines["right"].set_visible(False)
		ax4.spines["bottom"].set_position(("outward", 15))
		ax4.spines["left"].set_position(("outward", 15))
		ax4.spines["left"].set_visible(False)
		ax4.get_yaxis().set_visible(False)
		ax4.set_xlim([0, time_on[-1]/60])
		ax4.set_xticks(list(gen_start_time_on / 60) + [time_on[-1] / 60])
		ax4.set_xticklabels(np.arange(len(gen_start_time_on) + 1))

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
