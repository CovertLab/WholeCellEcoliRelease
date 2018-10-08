"""
Plots fold change for molecules related to tRNA charging (charged/uncharged tRNA
and synthetases) along with ribosome elongation rate.

Useful for seeing if any molecules have deviated far from initial values and
are possible causes for changes to the ribosomal elongation rate if the rate is
specified by tRNA charging in the simulation. tRNA and synthetases are shown
on a per amino acid basis with a sum from all relevant species.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/20/18

TODO: add amino acids and other metabolites involved in charging reactions
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.analysis.plotting_tools import COLORS_SMALL
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils.sparkline import whitePadSparklineAxis


def plot_ax(ax, x, y):
	'''
	Plots data and sets some common features for the axes

	Inputs:
		ax (matplotlib axes): axes to plot on
		x (numpy array of floats): x values to plot
		y (numpy array (1D or 2D) of floats): y values to plot
	'''

	ax.set_prop_cycle('color', COLORS_SMALL)  # Provides same color for each trace in different generations
	ax.plot(x, y)

def post_plot_formatting(ax, division_times, y_label, draw_horizontal=None, y_lim=None, show_x_axis=False):
	'''
	Formats axes after all data has been plotted

	Inputs:
		ax (matplotlib axes): axes to plot on
		division_times (numpy array of floats): times of division
		y_label (str): label for the y axis
		draw_horizontal (float): if not None, draws a horizontal line at the given y position
		y_lim (numpy array of floats): if not None, specifies the lower and upper y limits
		show_x_axis (bool): if True, displays x axis
	'''

	if y_lim is not None:
		ax.set_ylim(y_lim)
	ax.set_xlim([0, division_times[-1]])
	whitePadSparklineAxis(ax, xAxis=show_x_axis)
	if show_x_axis:
		ax.set_xticks(division_times[:-1], minor=True)
		ax.set_xticks([0, division_times[-1]])

	if draw_horizontal is not None:
		ax.axhline(draw_horizontal, color='k', linestyle='--', linewidth=1)

	str_format = FormatStrFormatter('%.3g')
	ax.xaxis.set_major_formatter(str_format)
	ax.yaxis.set_major_formatter(str_format)
	ax.set_ylabel(y_label, fontsize=8)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, 'seedOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		transcription = sim_data.process.transcription
		synthetase_names = transcription.synthetase_names
		uncharged_trna_names = transcription.rnaData['id'][transcription.rnaData['isTRna']]
		charged_trna_names = transcription.charged_trna_names
		aa_from_synthetase = transcription.aa_from_synthetase.T
		aa_from_trna = transcription.aa_from_trna.T

		aa_ids = sim_data.moleculeGroups.aaIDs
		n_aas = len(aa_ids)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)

		# Create plot and axes
		n_subplots = 6
		fig = plt.figure(figsize=(5, 10))
		rib_ax = plt.subplot(n_subplots, 1, 1)
		synth_ax = plt.subplot(n_subplots, 1, 2)
		frac_ax = plt.subplot(n_subplots, 1, 3)
		uncharged_trna_ax = plt.subplot(n_subplots, 1, 4)
		charged_trna_ax = plt.subplot(n_subplots, 1, 5)
		legend_ax = plt.subplot(n_subplots, 1, 6)

		initial_synthetase_counts = None
		initial_uncharged_trna_counts = None
		initial_charged_trna_counts = None
		division_times = []
		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))
			ribosome_reader = TableReader(os.path.join(simOutDir, "RibosomeData"))

			# Load data
			time = main_reader.readColumn('time') / 3600
			division_times.append(time[-1])
			elong_rate = ribosome_reader.readColumn("effectiveElongationRate")

			(synthetase_counts, uncharged_trna_counts, charged_trna_counts
				) = read_bulk_molecule_counts(simOutDir,
				(synthetase_names, uncharged_trna_names, charged_trna_names))

			## Synthetase counts
			synthetase_counts = np.dot(synthetase_counts, aa_from_synthetase)
			if initial_synthetase_counts is None:
				initial_synthetase_counts = synthetase_counts[1, :]
			normalized_synthetase_counts = synthetase_counts / initial_synthetase_counts

			## Uncharged tRNA counts
			uncharged_trna_counts = np.dot(uncharged_trna_counts, aa_from_trna)
			if initial_uncharged_trna_counts is None:
				initial_uncharged_trna_counts = uncharged_trna_counts[1, :]
			normalized_uncharged_trna_counts = uncharged_trna_counts / initial_uncharged_trna_counts

			## Charged tRNA counts
			charged_trna_counts = np.dot(charged_trna_counts, aa_from_trna)
			if initial_charged_trna_counts is None:
				initial_charged_trna_counts = charged_trna_counts[1, :]
			normalized_charged_trna_counts = charged_trna_counts / initial_charged_trna_counts

			## Fraction charged
			fraction_charged = charged_trna_counts / (charged_trna_counts + uncharged_trna_counts)

			# Plot data
			plot_ax(rib_ax, time[1:], elong_rate[1:])  # [1:] to remove spike down
			plot_ax(synth_ax, time, np.log2(normalized_synthetase_counts))
			plot_ax(frac_ax, time, fraction_charged)
			plot_ax(uncharged_trna_ax, time, np.log2(normalized_uncharged_trna_counts))
			plot_ax(charged_trna_ax, time, np.log2(normalized_charged_trna_counts))

		# Format plot axes
		post_plot_formatting(rib_ax, division_times, 'Ribosome\nElongation Rate', y_lim=[0, 22])
		post_plot_formatting(synth_ax, division_times, 'Synthetase\nFold Change', draw_horizontal=0)
		post_plot_formatting(frac_ax, division_times, 'Fraction\ntRNA Charged', y_lim=[0, 1])
		post_plot_formatting(uncharged_trna_ax, division_times, 'Uncharged tRNA\nFold Change', draw_horizontal=0)
		post_plot_formatting(charged_trna_ax, division_times, 'Charged tRNA\nFold Change', draw_horizontal=0, show_x_axis=True)
		charged_trna_ax.set_xlabel('Time (hr)')

		# Format and display legend below all plots
		legend_ax.set_prop_cycle('color', COLORS_SMALL)
		legend_ax.axis('off')
		legend_ax.plot(0, np.zeros((1, n_aas)))
		legend_ax.legend(aa_ids, ncol=3, fontsize=6, loc=10)

		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
