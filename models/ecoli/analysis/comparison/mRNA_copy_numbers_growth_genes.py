"""
Generates a comparison scatter plot of mRNA copy numbers for growth genes (genes
encoding for RNA polymerase or ribosomal subunits) from two sets of simulations,
one without operons and one with operons.
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


FIGSIZE = (4, 3.9)
BOUNDS = [0.5, 2.5]


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, validation_data2 = self.setup(input_sim_dir)

		cell_paths = ap2.get_cells()
		simOutDir = os.path.join(cell_paths[0], 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Get mask for mRNA genes that encode for ribosomal proteins or RNAPs
		all_cistron_ids = sim_data2.process.transcription.cistron_data['id']
		is_mRNA = sim_data2.process.transcription.cistron_data['is_mRNA']
		assert np.all(
			mRNA_cistron_ids == all_cistron_ids[is_mRNA])

		mRNA_is_rnap = sim_data2.process.transcription.cistron_data['is_RNAP'][is_mRNA]
		mRNA_is_ribosomal_protein = sim_data2.process.transcription.cistron_data['is_ribosomal_protein'][is_mRNA]
		all_plotted_mRNA_mask = np.logical_or(
			mRNA_is_rnap, mRNA_is_ribosomal_protein
			)

		def read_data(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			# Sample initial mRNA counts from each cell
			all_initial_counts = read_stacked_columns(
				cell_paths, 'mRNACounts', 'mRNA_cistron_counts',
				ignore_exception=True, fun=lambda x: x[0])

			return all_initial_counts

		c1 = read_data(ap1)
		c2 = read_data(ap2)

		if len(c1) == 0 or len(c2) == 0:
			print('Skipping analysis -- not enough sims run.')
			return

		m1 = c1.mean(axis=0)
		m2 = c2.mean(axis=0)

		# Normalize counts from two conditions
		ratio = m1[all_plotted_mRNA_mask].sum()/m2[all_plotted_mRNA_mask].sum()
		m2 = ratio * m2.astype(np.float)

		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)
		ax.plot(BOUNDS, BOUNDS, ls='--', lw=2, c='k', alpha=0.05)
		ax.scatter(
			np.log10(m1[mRNA_is_ribosomal_protein] + 1),
			np.log10(m2[mRNA_is_ribosomal_protein] + 1),
			c='turquoise', edgecolor='none', s=12, alpha=0.7,
			label=f'ribosomal subunits',
			clip_on=False)
		ax.scatter(
			np.log10(m1[mRNA_is_rnap] + 1),
			np.log10(m2[mRNA_is_rnap] + 1),
			c='darkslategray', edgecolor='none', s=12, alpha=0.7,
			label='RNAP subunits',
			clip_on=False)

		ax.set_xlabel('$\log_{10}$(mRNA copies + 1), old sims')
		ax.set_ylabel('$\log_{10}$(mRNA copies + 1), new sims')
		ax.set_xticks(np.arange(BOUNDS[0], BOUNDS[1] + 0.5, 0.5))
		ax.set_yticks(np.arange(BOUNDS[0], BOUNDS[1] + 0.5, 0.5))
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_position(("outward", 15))
		ax.spines["left"].set_position(("outward", 15))
		ax.set_xlim(BOUNDS)
		ax.set_ylim(BOUNDS)
		ax.legend(loc=2, prop={'size': 8})

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
