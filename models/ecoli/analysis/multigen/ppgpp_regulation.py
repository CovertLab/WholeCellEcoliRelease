"""
Analysis of ppGpp control and synthesis/degradation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/18/19
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import itertools
import os

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


RELA_RNA = 'EG10835_RNA[c]'
SPOT_RNA = 'EG10966_RNA[c]'


def read_data(cell_paths, table, column):
	sim_columns = []

	for sim_dir in cell_paths:
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		reader = TableReader(os.path.join(sim_out_dir, table))
		sim_columns.append(reader.readColumn2D(column))

	return np.vstack(sim_columns)

def remove_axes(ax, show_xaxis=False):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.tick_params(right=False, top=False)
	if not show_xaxis:
		ax.spines['bottom'].set_visible(False)
		ax.tick_params(bottom=False, labelbottom=False)

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# sim_data values
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		rna_data = sim_data.process.transcription.rna_data
		fractions = [
			'is_rRNA',
			'is_tRNA',
			'is_mRNA',
			'is_ribosomal_protein',
			'is_RNAP',
			]
		n_fractions = len(fractions)
		synthase_rna_idx_all_rnas = np.array([
			np.where(rna_data['id'] == RELA_RNA)[0][0],
			np.where(rna_data['id'] == SPOT_RNA)[0][0],
			])
		mrna_ids = rna_data['id'][rna_data['is_mRNA']]
		synthase_rna_idx_mrnas = np.array([
			np.where(mrna_ids == RELA_RNA)[0][0],
			np.where(mrna_ids == SPOT_RNA)[0][0],
			])
		synthase_order = ['relA', 'spoT']

		# Load simulation output
		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		cell_paths = ap.get_cells()
		time = read_data(cell_paths, 'Main', 'time').squeeze() / 60  # min
		counts_to_molar = read_data(cell_paths, 'EnzymeKinetics', 'countsToMolar').squeeze()
		synth_prob = read_data(cell_paths, 'RnaSynthProb', 'rnaSynthProb')
		reaction_rates = np.hstack((
			read_data(cell_paths, 'GrowthLimits', 'rela_syn'),
			read_data(cell_paths, 'GrowthLimits', 'spot_syn'),
			read_data(cell_paths, 'GrowthLimits', 'spot_deg'),
			))
		cell_count_iters = [
			read_bulk_molecule_counts(
				os.path.join(p, 'simOut'), [sim_data.molecule_ids.ppGpp])
			for p in cell_paths]
		ppgpp_count = np.hstack(list(itertools.chain.from_iterable(cell_count_iters)))
		del cell_count_iters
		ppgpp_conc = ppgpp_count * counts_to_molar * 1000  # uM
		mrna_count = read_data(cell_paths, 'mRNACounts', 'mRNA_counts')
		synthase_counts = mrna_count[:, synthase_rna_idx_mrnas]

		extra_plots = 4  # traces in addition to fractions
		grid_spec = GridSpec(n_fractions+extra_plots, 1)
		plt.figure(figsize=(5, 15))

		# ppGpp concentration
		mean_ppgpp = ppgpp_conc.mean()
		ax = plt.subplot(grid_spec[0, 0])
		ax.plot(time, ppgpp_conc)
		ax.axhline(mean_ppgpp, linestyle='--', color='k')
		ax.set_ylabel('ppGpp conc (uM)', fontsize=10)
		remove_axes(ax)

		# ppGpp reaction rates
		total_rate = reaction_rates[:, 0] + reaction_rates[:, 1] - reaction_rates[:, 2]
		ax = plt.subplot(grid_spec[1, 0])
		ax.plot(time, reaction_rates, alpha=0.8)
		ax.plot(time, total_rate, 'k', alpha=0.3, linewidth=1.)
		ax.legend(['RelA syn', 'SpoT syn', 'SpoT deg', 'Total'], fontsize=8)
		ax.set_ylabel('ppGpp reaction rates\n(uM/s)', fontsize=10)
		remove_axes(ax)

		# Synthase probabilities
		ax = plt.subplot(grid_spec[2, 0])
		ax.plot(time, synth_prob[:, synthase_rna_idx_all_rnas], alpha=0.8)
		ax.legend(synthase_order, fontsize=8)
		ax.set_ylabel('RNA synth prob', fontsize=10)
		remove_axes(ax)

		# Synthase RNA counts
		ax = plt.subplot(grid_spec[3, 0])
		ax.plot(time, synthase_counts, alpha=0.8)
		ax.legend(synthase_order, fontsize=8)
		ax.set_ylabel('RNA counts', fontsize=10)
		remove_axes(ax)

		# Sum of probabilities for RNA fractions
		for i, fraction in enumerate(fractions):
			ax = plt.subplot(grid_spec[i+extra_plots, 0])
			ax.plot(time, synth_prob[:, rna_data[fraction]].sum(1))
			ax.set_ylabel(fraction, fontsize=10)
			remove_axes(ax, show_xaxis=(i == n_fractions - 1))

		# TODO: add all positive/negative regulated expression

		ax.set_xlabel('Time (min)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
