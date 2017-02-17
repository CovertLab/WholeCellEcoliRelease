#!/usr/bin/env python
"""
Compare fluxes in simulation to target fluxes

@date: Created 12/15/16
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os
import cPickle
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.plotting_tools import COLORS_LARGE

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

BURN_IN_STEPS = 20

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile))

	constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly
	constrainedReactions = np.array(sim_data.process.metabolism.constrainedReactionList)

	mainListener = TableReader(os.path.join(simOutDir, "Main"))
	initialTime = mainListener.readAttribute("initialTime")
	time = mainListener.readColumn("time") - initialTime
	timeStepSec = mainListener.readColumn("timeStepSec")
	mainListener.close()

	massListener = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = massListener.readColumn("cellMass")
	dryMass = massListener.readColumn("dryMass")
	massListener.close()

	coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(units.g / units.L) * timeStepSec # units - g.s/L

	# read constraint data
	enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	targetFluxes = (units.dmol / units.g / units.s) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
	actualFluxes = (units.dmol / units.g / units.s) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
	reactionConstraint = enzymeKineticsReader.readColumn("reactionConstraint")
	enzymeKineticsReader.close()

	targetFluxes = targetFluxes.asNumber(units.mmol / units.g / units.h)
	actualFluxes = actualFluxes.asNumber(units.mmol / units.g / units.h)

	targetAve = np.mean(targetFluxes[BURN_IN_STEPS:, :], axis = 0)
	actualAve = np.mean(actualFluxes[BURN_IN_STEPS:, :], axis = 0)

	relError = np.abs((actualFluxes[BURN_IN_STEPS:, :] - targetFluxes[BURN_IN_STEPS:, :]) / (targetFluxes[BURN_IN_STEPS:, :] + 1e-15))
	aveError = np.mean(relError, axis = 0)

	kcatOnlyReactions = np.all(constraintIsKcatOnly[reactionConstraint[BURN_IN_STEPS:,:]], axis = 0)
	kmAndKcatReactions = ~np.any(constraintIsKcatOnly[reactionConstraint[BURN_IN_STEPS:,:]], axis = 0)
	mixedReactions = ~(kcatOnlyReactions ^ kmAndKcatReactions)

	kmAndKcatThresholds = [2, 10]
	kmAndKcatCategorization = np.zeros(np.sum(kmAndKcatReactions))
	categorization = np.zeros(reactionConstraint.shape[1])
	for i, threshold in enumerate(kmAndKcatThresholds):
		# kmAndKcatCategorization[targetAve[kmAndKcatReactions] / actualAve[kmAndKcatReactions] > threshold] = i + 1
		kmAndKcatCategorization[actualAve[kmAndKcatReactions] / targetAve[kmAndKcatReactions] > threshold] = i + 1
		categorization[actualAve / targetAve > threshold] = i + 1
	kmAndKcatCategorization[actualAve[kmAndKcatReactions] == 0] = -1

	kcatOnlyThresholds = [2, 10]
	kcatOnlyCategorization = np.zeros(np.sum(kcatOnlyReactions))
	for i, threshold in enumerate(kcatOnlyThresholds):
		kcatOnlyCategorization[actualAve[kcatOnlyReactions] / targetAve[kcatOnlyReactions] > threshold] = i + 1
	kcatOnlyCategorization[actualAve[kcatOnlyReactions] == 0] = -1

	# url for ecocyc to highlight fluxes that are 0 on metabolic network diagram
	siteStr = "https://ecocyc.org/overviewsWeb/celOv.shtml?zoomlevel=1&orgid=ECOLI"
	excluded = ['RXN0-2201', 'RXN-16000', 'RXN-12583', 'RXN-11496', 'DIMESULFREDUCT-RXN', '3.6.1.41-R[4/63051]5-NUCLEOTID-RXN'] # reactions not recognized by ecocyc
	rxns = []
	for i, reaction in enumerate(constrainedReactions):
		if actualAve[i] == 0:
			rxn = re.findall(".+RXN", reaction)
			if len(rxn) == 0:
				rxn = re.findall("RXN[^-]*-[0-9]+", reaction)
			if rxn[0] not in excluded:
				siteStr += "&rnids=%s" % rxn[0]
			rxns.append(rxn[0])
	# print siteStr

	csvFile = open(os.path.join(plotOutDir, plotOutFileName + ".tsv"), "wb")
	output = csv.writer(csvFile, delimiter = "\t")
	output.writerow(["ecocyc link:", siteStr])
	output.writerow(["kM and kcat", "Target", "Actual", "Category"])
	for reaction, target, flux, category in zip(constrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], kmAndKcatCategorization):
		output.writerow([reaction, target, flux, category])

	output.writerow(["kcat only"])
	for reaction, target, flux, category in zip(constrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], kcatOnlyCategorization):
		output.writerow([reaction, target, flux, category])

	csvFile.close()

	targetAve[targetAve == 0] += 1e-9
	actualAve[actualAve == 0] += 1e-9

	plt.figure(figsize = (8, 8))
	from scipy.stats import pearsonr
	targetPearson = targetAve[kmAndKcatReactions]
	actualPearson = actualAve[kmAndKcatReactions]
	# plt.title(pearsonr(np.log10(targetPearson[actualPearson > 0]), np.log10(actualPearson[actualPearson > 0])))
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == 0], actualAve[kmAndKcatReactions][kmAndKcatCategorization == 0], "og")
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == 1], actualAve[kmAndKcatReactions][kmAndKcatCategorization == 1], "o")
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == 2], actualAve[kmAndKcatReactions][kmAndKcatCategorization == 2], "or")
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == -1], actualAve[kmAndKcatReactions][kmAndKcatCategorization == -1], "og")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 0], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 0], "og")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 1], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 1], "o")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 2], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 2], "or")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == -1], actualAve[kcatOnlyReactions][kcatOnlyCategorization == -1], "og")
	plt.loglog(targetAve[categorization == 0], actualAve[categorization == 0], "og")
	plt.loglog(targetAve[categorization == 1], actualAve[categorization == 1], "o")
	plt.loglog(targetAve[categorization == 2], actualAve[categorization == 2], "or")
	plt.loglog(targetAve[categorization == -1], actualAve[categorization == -1], "og")
	# plt.loglog(targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], "o")
	# plt.loglog(targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], "ro")
	plt.loglog([1e-10, 1e4], [1e-10, 1e4], '--g')
	plt.loglog([1e-10, 1e3], [1e-9, 1e4], '--r')
	# plt.loglog([1e-13, 1], [1e-14, 0.1], '--r')
	plt.xlabel("Target Flux (mmol/g/hr)")
	plt.ylabel("Actual Flux (mmol/g/hr)")
	plt.minorticks_off()
	whitePadSparklineAxis(plt.axes())

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

	# Bokeh
	from bokeh.plotting import figure, output_file, ColumnDataSource, show
	from bokeh.charts import Bar
	from bokeh.models import HoverTool, BoxZoomTool, LassoSelectTool, PanTool, WheelZoomTool, ResizeTool, UndoTool, RedoTool

	source = ColumnDataSource(
		data = dict(
			x = targetAve,
			y = actualAve,
			reactionName = constrainedReactions)
		)

	hover = HoverTool(
		tooltips = [
			("Reaction", "@reactionName"),
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

	p1 = figure(x_axis_label = "Target",
		x_axis_type = "log",
		x_range = [min(targetAve[targetAve > 0]), max(targetAve)],
		y_axis_label = "Actual",
		y_axis_type = "log",
		y_range = [min(actualAve[actualAve > 0]), max(actualAve)],
		width = 800,
		height = 800,
		tools = TOOLS,
		)

	p1.scatter(targetAve, actualAve, source = source, size = 8)
	p1.line([1e-15, 10], [1e-15, 10], line_color = "red", line_dash = "dashed")
	

	## bar plot of error
	# sortedReactions = [constrainedReactions[x] for x in np.argsort(aveError)[::-1]]
	# aveError[np.log10(aveError) == -np.inf] = 0

	# source = ColumnDataSource(
	# 	data = dict(
	# 			x = sorted(relError, reverse = True),
	# 			reactionName = sortedReactions
	# 		)
	# 	)

	# p2 = Bar(data, values = "x")

	# hover2 = p2.select(dict(type=HoverTool))
	# hover2.tooltips = [("Reaction", "@reactionName")]

	## flux for each reaction
	hover2 = HoverTool(
		tooltips = [
			("Reaction", "@reactionName"),
			]
		)

	TOOLS2 = [hover2, 
		BoxZoomTool(),
		LassoSelectTool(),
		PanTool(),
		WheelZoomTool(),
		ResizeTool(),
		UndoTool(),
		RedoTool(),
		 "reset"
		 ]

	p2 = figure(x_axis_label = "Time(s)", 
		y_axis_label = "Flux",
		y_axis_type = "log",
		y_range = [1e-8, 1],
		width = 800,
		height = 800,
		tools = TOOLS2,
		)

	colors = COLORS_LARGE
	nTimesteps = len(time[BURN_IN_STEPS:])
	x = time[BURN_IN_STEPS:]
	y = actualFluxes[BURN_IN_STEPS:, 0]
	reactionName = np.repeat(constrainedReactions[0], nTimesteps)

	source = ColumnDataSource(
		data = dict(
			x = x,
			y = y,
			reactionName = reactionName)
		)

	p2.line(x, y, line_color = colors[0], source = source)

	# Plot remaining metabolites onto initialized figure
	for m in np.arange(1, actualFluxes.shape[1]):
		y = actualFluxes[BURN_IN_STEPS:, m]
		reactionName = np.repeat(constrainedReactions[m], nTimesteps)

		source = ColumnDataSource(
			data = dict(
				x = x,
				y = y,
				reactionName = reactionName)
		)

		p2.line(x, y, line_color = colors[m % len(colors)], source = source)

	if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
		os.makedirs(os.path.join(plotOutDir, "html_plots"))

	import bokeh.io
	p = bokeh.io.vplot(p1, p2)
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
