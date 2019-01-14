"""
Shows fold change of metabolites over the course of the simulation

@date: Created 8/25/2016
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import bokeh.io
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool

from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

PLOT_BOKEH = False


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		metaboliteNames = enzymeKineticsdata.readAttribute("metaboliteNames")
		metaboliteCounts = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
		normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]

		# Read time info from the listener
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		enzymeKineticsdata.close()

		colors = COLORS_LARGE # to match colors between the pdf and html plots
		plt.figure(figsize = (8.5, 11))
		ax = plt.subplot(1, 1, 1)
		ax.set_prop_cycle('color', colors)
		plt.plot(time, normalizedCounts)
		plt.xlabel("Time (s)")
		plt.ylabel("Metabolite Fold Change")
		# plt.legend(metaboliteNames[lowCountsIdx], fontsize = 4)
		# plt.plot(time, normalizedCounts[:, lowCountsIdx])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		if not PLOT_BOKEH:
			return

		if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
			os.makedirs(os.path.join(plotOutDir, "html_plots"))
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
