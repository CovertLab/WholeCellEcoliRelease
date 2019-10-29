"""
Assess the sensitivity of succinate and isocitrate dehydrogenase fluxes to
each kinetic constraint being disabled.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/25/19
"""

from __future__ import absolute_import
from __future__ import division

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
import numpy as np
import os

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath, sparkline


def calc_z(data):
	"""
	Calculates a z score for relative flux time series data.

	Args:
		data (ndarray[float]): 2D array (n time steps x m reactions) of fluxes
			with the last column being the reference flux

	Returns:
		ndarray[float]: z score for each reaction
	"""

	rel_flux = (data[:, :-1] - data[:, -1:]) / data[:, -1:]
	rel_mean = np.sort(rel_flux, axis=0)[:5, :].mean(axis=0)
	return (rel_mean - rel_mean.mean()) / rel_mean.std()

def plot_lows(ax, data, threshold, label):
	"""
	Plots the largest negative change for each reaction.

	Args:
		ax (matplotlib.axes): axes to plot on
		data (ndarray[float]): 2D array (n time steps x m reactions) of fluxes
			with the last column being the reference flux
		threshold (float): threshold to plot
		label (str): reaction label for y axis
	"""

	# Plot data
	ax.set_yscale('symlog', threshold=0.01)
	ax.fill_between(range(len(data)), sorted(data), color='b')
	ax.axhline(threshold, color='k', linestyle='--', linewidth=0.5)

	# Format axes
	y_ticks = np.hstack((ax.get_yticks(), threshold))
	sparkline.whitePadSparklineAxis(ax, xAxis=False)
	ax.set_xticks([])
	ax.set_yticks(y_ticks)
	ax.tick_params(labelsize=6)
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax.set_ylabel('Maximum change in\n{} flux'.format(label), fontsize=6)

def plot_threshold(ax, data, threshold, reactions):
	"""
	Plots the largest negative change for each reaction the changes the flux
	more than the threshold.  Prints the full reaction names for these reactions.

	Args:
		ax (matplotlib.axes): axes to plot on
		data (ndarray[float]): 2D array (n time steps x m reactions) of fluxes
			with the last column being the reference flux
		threshold (float): threshold to plot and use as cutoff
		reactions (ndarray[str]): names for each reaction in data
	"""

	sorted_idx = np.argsort(data)
	below_idx = np.where(data[sorted_idx] < threshold)[0]

	# Plot data
	ax.set_yscale('symlog', threshold=0.01)
	ax.bar(below_idx, data[sorted_idx[below_idx]], color='b')
	ax.axhline(threshold, color='k', linestyle='--', linewidth=0.5)

	# Format axes
	y_ticks = np.hstack((ax.get_yticks(), threshold))
	sparkline.whitePadSparklineAxis(ax)
	ax.spines["bottom"].set_visible(False)
	ax.tick_params(bottom=False)
	ax.set_xticks(below_idx)
	ax.set_xticklabels([r[:10] for r in reactions[sorted_idx[below_idx]]], rotation=90)
	ax.set_yticks(y_ticks)
	ax.tick_params(labelsize=6)
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

	print(reactions[sorted_idx[below_idx]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata.get('variant', '') != 'flux_sensitivity':
			print 'This plot only runs for the flux_sensitivity variant.'
			return

		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		succ_fluxes = []
		iso_fluxes = []
		for variant in variants:
			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				fba_reader = TableReader(os.path.join(simOutDir, 'FBAResults'))

				# Load data
				reactions = np.array(fba_reader.readAttribute('sensitivity_reactions'))
				succ_fluxes += [fba_reader.readColumn('succinate_flux_sensitivity')[1:, :]]
				iso_fluxes += [fba_reader.readColumn('isocitrate_flux_sensitivity')[1:, :]]

		succ_fluxes = np.vstack(succ_fluxes)
		iso_fluxes = np.vstack(iso_fluxes)

		succ_z = calc_z(succ_fluxes)
		iso_z = calc_z(iso_fluxes)

		threshold = -0.1

		# Plot data
		plt.figure()
		gs = gridspec.GridSpec(2, 2)

		## Succinate dehydrogenase all fluxes
		ax = plt.subplot(gs[0, 0])
		plot_lows(ax, succ_z, threshold, 'succinate dehydrogenase')

		## Succinate dehydrogenase fluxes over threshold
		ax = plt.subplot(gs[0, 1])
		plot_threshold(ax, succ_z, threshold, reactions)

		## Isocitrate dehydrogenase all fluxes
		ax = plt.subplot(gs[1, 0])
		plot_lows(ax, iso_z, threshold, 'isocitrate dehydrogenase')

		## Isocitrate dehydrogenase fluxes over threshold
		ax = plt.subplot(gs[1, 1])
		plot_threshold(ax, iso_z, threshold, reactions)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
