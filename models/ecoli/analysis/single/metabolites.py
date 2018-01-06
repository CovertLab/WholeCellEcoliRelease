#!/usr/bin/env python
"""
Shows fold change of metabolites over the course of the simulation

@date: Created 8/25/2016
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from random import random

from wholecell.utils import units
import cPickle
import ast
import itertools

from wholecell.utils.modular_fba import FluxBalanceAnalysis

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.analysis.plotting_tools import COLORS_LARGE, COLORS_SMALL

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

PLOT_BOKEH = False

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	metaboliteNames = enzymeKineticsdata.readAttribute("metaboliteNames")
	metaboliteCounts = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
	normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
	
	# Read time info from the listener
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	enzymeKineticsdata.close()

	colors = COLORS_LARGE # to match colors between the pdf and html plots
	plt.figure(figsize = (8.5, 11))
	ax = plt.subplot(1, 1, 1)
	ax.set_color_cycle(colors)
	plt.plot(time, normalizedCounts)
	plt.xlabel("Time (s)")
	plt.ylabel("Metabolite Fold Change")
	# plt.legend(metaboliteNames[lowCountsIdx], fontsize = 4)
	# plt.plot(time, normalizedCounts[:, lowCountsIdx])

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

	if not PLOT_BOKEH:
		return

	# Bokeh
	from bokeh.plotting import figure, output_file, ColumnDataSource, show
	from bokeh.models import HoverTool, BoxZoomTool, LassoSelectTool, PanTool, WheelZoomTool, ResizeTool, UndoTool, RedoTool
	import bokeh.io
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

	import bokeh.io
	bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + ".html"), title = plotOutFileName, autosave = False)
	bokeh.io.save(plot)
	bokeh.io.curstate().reset()

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
