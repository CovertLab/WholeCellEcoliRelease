"""
Compare metabolite concentrations to target concentrations in FBA objective

@date: Created 2/6/18
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


def plot_data(gs, col, time, x, y, label, molecule_names):
	x_ave = np.nanmean(x, axis=0)
	y_ave = np.nanmean(y, axis=0)

	# Plot average comparison with lines denoting order of magnitude
	ax = plt.subplot(gs[0, col])
	ax.plot([-6, 6], [-6, 6], 'k')
	ax.plot([-5, 6], [-6, 5], 'k')
	ax.plot([-6, 5], [-5, 6], 'k')
	ax.plot(np.log10(x_ave), np.log10(y_ave), "ob", markeredgewidth=0,
	        alpha=0.25)
	ax.set_xlabel("Log10(Target Concentration [mmol/L])", fontsize=8)
	ax.set_ylabel("Log10({} Concentration [mmol/L])".format(label), fontsize=8)
	ax.set_title(label)

	whitePadSparklineAxis(ax)
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	ax.set_ylim(ylim[0] - 0.1, ylim[1])
	ax.set_xlim(xlim[0] - 0.1, xlim[1])
	ax.set_yticks(range(-6, int(ylim[1]) + 1, 2))
	ax.set_xticks(range(-6, int(xlim[1]) + 1, 2))
	ax.tick_params(axis='both', which='major', labelsize=6)

	# Plot ratio of actual concentration to target concentration
	ratio = np.log10(y / x)
	ax = plt.subplot(gs[1, col])
	ax.plot(time / 60, ratio)
	ax.set_xlabel("Time (min)", fontsize=8)
	ax.set_ylabel("Log10(Concentration to Target)", fontsize=8)
	ax.tick_params(axis='both', which='major', labelsize=6)

	# Plot outliers of ratio
	means = np.mean(ratio, axis=0)
	mean = np.mean(means[np.isfinite(means)])
	std = np.std(means[np.isfinite(means)])
	outliers = np.unique(
		np.where((ratio[1:, :] > mean + std / 2) | (ratio[1:, :] < mean - std / 2))[1])
	outliers = outliers[np.argsort(np.abs(means[outliers]))][::-1][:10]  # limit outliers to 10 with largest mean
	idx = outliers[np.argsort(means[outliers])][::-1]

	ax = plt.subplot(gs[2, col])
	if len(idx):
		ax.plot(time / 60, ratio[:, idx])
		ax.legend(molecule_names[idx], fontsize=6)
	ax.set_xlabel("Time (min)", fontsize=8)
	ax.set_ylabel("Log10(Concentration to Target) Outliers", fontsize=8)
	ax.tick_params(axis='both', which='major', labelsize=6)


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Read data from listeners
		enzymeKinetics = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		metabolite_names = enzymeKinetics.readAttribute('metaboliteNames')
		initial_counts = enzymeKinetics.readColumn("metaboliteCountsInit")[1:, :]
		counts_to_molar = enzymeKinetics.readColumn('countsToMolar')[1:].reshape(-1, 1)
		initial_conc = counts_to_molar * initial_counts

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		target_conc = fbaResults.readColumn("targetConcentrations")[1:,:]
		molecule_names = np.array(fbaResults.readAttribute("homeostaticTargetMolecules"))

		mainReader = TableReader(os.path.join(simOutDir, "Main"))
		time = mainReader.readColumn("time")[1:] - mainReader.readAttribute("initialTime")

		actual_counts, = read_bulk_molecule_counts(simOutDir, (metabolite_names,))
		actual_conc = actual_counts[1:, :] * counts_to_molar

		plt.figure(figsize = (8.5, 11))
		gs = gridspec.GridSpec(3, 2)
		plot_data(gs, 0, time, target_conc, initial_conc, 'Initial', molecule_names)
		plot_data(gs, 1, time, target_conc, actual_conc, 'Actual', molecule_names)

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
