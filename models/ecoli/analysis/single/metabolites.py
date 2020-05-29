"""
Shows fold change of metabolites over the course of the simulation

@date: Created 8/25/2016
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

import numpy as np
from matplotlib import pyplot as plt
import bokeh.io
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool

from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import filepath
from models.ecoli.analysis import singleAnalysisPlot


PLOT_BOKEH = False


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile) as f:
			sim_data = cPickle.load(f)
		aa_ids = sim_data.moleculeGroups.aaIDs

		# Listeners used
		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		main_reader = TableReader(os.path.join(simOutDir, "Main"))

		# Metabolite data
		metaboliteNames = np.array(enzymeKineticsdata.readAttribute("metaboliteNames"))
		metaboliteCounts = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
		normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
		aa_mask = np.array([m in aa_ids for m in metaboliteNames])

		# Highlight outliers
		mean_final = normalizedCounts[-1, ~aa_mask].mean()
		std_final = normalizedCounts[-1, ~aa_mask].std()
		highlighted = (
			(normalizedCounts[-1, :] > mean_final + 3*std_final)
			| (normalizedCounts[-1, :] < mean_final - 3*std_final)
		)

		# Sort amino acids for labeling
		sorted_idx = np.argsort(normalizedCounts[-1, aa_mask])[::-1]

		# Read time info from the listener
		initialTime = main_reader.readAttribute("initialTime")
		time = (main_reader.readColumn("time") - initialTime) / 60

		colors = COLORS_LARGE
		plt.figure(figsize = (8.5, 11))

		# Plot everything but amino acids
		ax = plt.subplot(2, 1, 1)
		ax.set_prop_cycle('color', colors)

		## Plot and label metabolites that are different from the mean
		mask = ~aa_mask & highlighted
		if np.any(mask):
			plt.plot(time, normalizedCounts[:, mask])
			plt.legend(metaboliteNames[mask], fontsize=8)

		## Plot the rest of the metabolites that are not amino acids
		mask = ~aa_mask & ~highlighted
		if np.any(mask):
			plt.plot(time, normalizedCounts[:, mask])

		## Formatting
		plt.xlabel("Time (min)")
		plt.ylabel("Metabolite fold change")
		plt.title('All metabolites (excluding amino acids)')

		# Plot only amino acids
		ax = plt.subplot(2, 1, 2)
		ax.set_prop_cycle('color', colors)
		plt.plot(time, normalizedCounts[:, aa_mask][:, sorted_idx])
		plt.legend(metaboliteNames[aa_mask][sorted_idx], fontsize=8, ncol=2)
		plt.xlabel("Time (min)")
		plt.ylabel("Metabolite fold change - amino acids only")
		plt.title('Only amino acids')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		if not PLOT_BOKEH:
			return

		filepath.makedirs(plotOutDir, "html_plots")
		bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + ".html"), title = plotOutFileName, autosave = False)

		nTimesteps = time.shape[0]
		nMolecules = normalizedCounts.shape[1]
		freq = 100
		x_a = time[::freq]

		plot = figure(x_axis_label = "Time(s)", y_axis_label = "Metabolite Fold Change", width = 800, height = 800)
		for m in np.arange(nMolecules):
			y = normalizedCounts[:, m]
			y_a = y[::freq]
			metaboliteName = np.repeat(metaboliteNames[m], x_a.shape)
			source = ColumnDataSource(data = dict(x = x_a, y = y_a, ID = metaboliteName))
			circle = plot.circle("x", "y", alpha = 0, source = source)
			plot.add_tools(HoverTool(tooltips = [("ID", "@ID")], renderers = [circle]))
			plot.line(time, y, line_color = colors[m % len(colors)])

		bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + ".html"), title = plotOutFileName, autosave = False)
		bokeh.io.save(plot)
		bokeh.io.curstate().reset()


if __name__ == "__main__":
	Plot().cli()
