"""
Analysis of ppGpp control and synthesis/degradation.
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_bulk_molecules, read_stacked_columns


RELA_RNA = 'EG10835_RNA'
SPOT_RNA = 'EG10966_RNA'


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# sim_data values
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		cistron_data = sim_data.process.transcription.cistron_data
		fractions = [
			'is_rRNA',
			'is_tRNA',
			'is_mRNA',
			'is_ribosomal_protein',
			'is_RNAP',
			]
		n_fractions = len(fractions)
		synthase_cistron_idx_all_cistrons = np.array([
			np.where(cistron_data['id'] == RELA_RNA)[0][0],
			np.where(cistron_data['id'] == SPOT_RNA)[0][0],
			])
		mrna_cistron_ids = cistron_data['id'][cistron_data['is_mRNA']]
		synthase_cistron_idx_mrnas = np.array([
			np.where(mrna_cistron_ids == RELA_RNA)[0][0],
			np.where(mrna_cistron_ids == SPOT_RNA)[0][0],
			])
		synthase_order = ['relA', 'spoT']

		# Load simulation output
		cell_paths = self.ap.get_cells()
		time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze() / 60  # min
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar').squeeze()
		synth_prob_per_cistron = read_stacked_columns(cell_paths, 'RnaSynthProb', 'rna_synth_prob_per_cistron')
		reaction_rates = np.vstack((
			read_stacked_columns(cell_paths, 'GrowthLimits', 'rela_syn').sum(axis=1),
			read_stacked_columns(cell_paths, 'GrowthLimits', 'spot_syn').squeeze(),
			read_stacked_columns(cell_paths, 'GrowthLimits', 'spot_deg').squeeze(),
			)).T
		ppgpp_count, = read_stacked_bulk_molecules(cell_paths, [sim_data.molecule_ids.ppGpp])
		ppgpp_conc = ppgpp_count * counts_to_molar * 1000  # uM
		mrna_cistron_count = read_stacked_columns(cell_paths, 'mRNACounts', 'mRNA_cistron_counts')
		synthase_counts = mrna_cistron_count[:, synthase_cistron_idx_mrnas]

		extra_plots = 4  # traces in addition to fractions
		grid_spec = GridSpec(n_fractions+extra_plots, 1)
		plt.figure(figsize=(5, 15))

		# ppGpp concentration
		mean_ppgpp = ppgpp_conc.mean()
		ax = plt.subplot(grid_spec[0, 0])
		ax.plot(time, ppgpp_conc)
		ax.axhline(mean_ppgpp, linestyle='--', color='k')
		ax.set_ylabel('ppGpp conc (uM)', fontsize=10)
		self.remove_border(ax=ax, bottom=True)

		# ppGpp reaction rates
		total_rate = reaction_rates[:, 0] + reaction_rates[:, 1] - reaction_rates[:, 2]
		ax = plt.subplot(grid_spec[1, 0])
		ax.plot(time, reaction_rates, alpha=0.8)
		ax.plot(time, total_rate, 'k', alpha=0.3, linewidth=1.)
		ax.legend(['RelA syn', 'SpoT syn', 'SpoT deg', 'Total'], fontsize=8)
		ax.set_ylabel('ppGpp reaction rates\n(uM/s)', fontsize=10)
		self.remove_border(ax=ax, bottom=True)

		# Synthase probabilities
		ax = plt.subplot(grid_spec[2, 0])
		ax.plot(time, synth_prob_per_cistron[:, synthase_cistron_idx_all_cistrons], alpha=0.8)
		ax.legend(synthase_order, fontsize=8)
		ax.set_ylabel('RNA synth prob', fontsize=10)
		self.remove_border(ax=ax, bottom=True)

		# Synthase RNA counts
		ax = plt.subplot(grid_spec[3, 0])
		ax.plot(time, synthase_counts, alpha=0.8)
		ax.legend(synthase_order, fontsize=8)
		ax.set_ylabel('RNA counts', fontsize=10)
		self.remove_border(ax=ax, bottom=True)

		# Sum of probabilities for RNA fractions
		for i, fraction in enumerate(fractions):
			ax = plt.subplot(grid_spec[i+extra_plots, 0])
			ax.plot(time, synth_prob_per_cistron[:, cistron_data[fraction]].sum(1))
			ax.set_ylabel(fraction, fontsize=10)
			self.remove_border(ax=ax, bottom=(i != n_fractions - 1))

		# TODO: add all positive/negative regulated expression

		ax.set_xlabel('Time (min)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
