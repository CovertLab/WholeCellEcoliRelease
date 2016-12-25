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

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

	metaboliteCounts = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
	normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
	
	# Read time info from the listener
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))

	# Load constants
	nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / COUNTS_UNITS)
	cellDensity = sim_data.constants.cellDensity.asNumber(MASS_UNITS/VOLUME_UNITS)

	objective = dict(
			(key, sim_data.process.metabolism.concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in sim_data.process.metabolism.concDict
			)

	externalExchangeMolecules = sim_data.nutrientData["secretionExchangeMolecules"]
	for t, nutrientsLabel in sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel]:
		externalExchangeMolecules += sim_data.nutrientData["importExchangeMolecules"][nutrientsLabel]
	externalExchangeMolecules = sorted(externalExchangeMolecules)
	extMoleculeMasses = sim_data.getter.getMass(externalExchangeMolecules)

	moleculeMasses = dict(zip(
		externalExchangeMolecules,
		sim_data.getter.getMass(externalExchangeMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)
		))

	fba = FluxBalanceAnalysis(
			sim_data.process.metabolism.reactionStoich,
			externalExchangeMolecules,
			objective,
			objectiveType = "standard",
			moleculeMasses = moleculeMasses,
			solver = 'glpk'
			)

	metaboliteNames = np.array(fba.outputMoleculeIDs())

	# get only fluxes that drop
	# lowCountsIdx = np.unique(np.where(normalizedCounts[1:, :] < 0.99)[1])
	# get only fluxes that don't double by end
	# lowCountsIdx = np.where(normalizedCounts[-1, :] < 1.5)[0]
	# get only fluxes that more than double by end
	# lowCountsIdx = np.where(normalizedCounts[-1, :] > 2.1)[0]

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


	# Bokeh
	from bokeh.plotting import figure, output_file, ColumnDataSource, show
	from bokeh.models import HoverTool, BoxZoomTool, LassoSelectTool, PanTool, WheelZoomTool, ResizeTool, UndoTool, RedoTool

	nTimesteps = time.shape[0]
	nMolecules = normalizedCounts.shape[1]

	# Only create Bokeh plot if every metabolite is identified
	if metaboliteNames.shape[0] != nMolecules:
		return

	moleculeCounts = np.round(10* np.random.rand(nTimesteps*nMolecules).reshape((nMolecules, nTimesteps)))

	# Plot first metabolite to initialize plot settings
	x = time
	y = normalizedCounts[:, 0]
	metaboliteName = np.repeat(metaboliteNames[0], nTimesteps)

	source = ColumnDataSource(
		data = dict(
			x = x,
			y = y,
			metaboliteName = metaboliteName)
		)

	hover = HoverTool(
		tooltips = [
			("ID", "@metaboliteName"),
			]
		)
	
	TOOLS = [hover, 
		BoxZoomTool(),
		LassoSelectTool(),
		PanTool(),
		WheelZoomTool(),
		ResizeTool(),
		UndoTool(),
		RedoTool(),
		 "reset"
		 ]

	p = figure(x_axis_label = "Time(s)", 
		y_axis_label = "Metabolite Fold Change",
		width = 800,
		height = 800,
		tools = TOOLS,
		)

	p.line(x, y, line_color = colors[0], source = source)

	# Plot remaining metabolites onto initialized figure
	for m in np.arange(1, nMolecules):
		x = time
		y = normalizedCounts[:, m]
		metaboliteName = np.repeat(metaboliteNames[m], nTimesteps)

		source = ColumnDataSource(
			data = dict(
				x = x,
				y = y,
				metaboliteName = metaboliteName))

		p.line(x, y, line_color = colors[m % len(colors)], source = source)

	if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
		os.makedirs(os.path.join(plotOutDir, "html_plots"))

	import bokeh.io
	bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + ".html"), title=plotOutFileName, autosave=False)
	bokeh.io.save(p)


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

	args = parser.parse_args().__dict__
	
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
