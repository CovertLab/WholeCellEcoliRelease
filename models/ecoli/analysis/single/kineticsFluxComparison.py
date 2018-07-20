"""
Compare fluxes in simulation to target fluxes

@date: Created 12/15/16
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle
import csv
import re

import numpy as np
from matplotlib import pyplot as plt
import bokeh.io
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import (HoverTool, BoxZoomTool, LassoSelectTool, PanTool,
	WheelZoomTool, ResizeTool, UndoTool, RedoTool)

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.plotting_tools import COLORS_LARGE

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

BURN_IN_STEPS = 20


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile))

		constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly

		mainListener = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = mainListener.readAttribute("initialTime")
		time = mainListener.readColumn("time") - initialTime
		mainListener.close()

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass")
		dryMass = massListener.readColumn("dryMass")
		massListener.close()

		coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		# read constraint data
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		targetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
		actualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
		reactionConstraint = enzymeKineticsReader.readColumn("reactionConstraint")
		constrainedReactions = np.array(enzymeKineticsReader.readAttribute("constrainedReactions"))
		enzymeKineticsReader.close()

		targetFluxes = targetFluxes.asNumber(units.mmol / units.g / units.h)
		actualFluxes = actualFluxes.asNumber(units.mmol / units.g / units.h)

		targetAve = np.mean(targetFluxes[BURN_IN_STEPS:, :], axis = 0)
		actualAve = np.mean(actualFluxes[BURN_IN_STEPS:, :], axis = 0)

		kcatOnlyReactions = np.all(constraintIsKcatOnly[reactionConstraint[BURN_IN_STEPS:,:]], axis = 0)
		kmAndKcatReactions = ~np.any(constraintIsKcatOnly[reactionConstraint[BURN_IN_STEPS:,:]], axis = 0)
		mixedReactions = ~(kcatOnlyReactions ^ kmAndKcatReactions)

		thresholds = [2, 10]
		categorization = np.zeros(reactionConstraint.shape[1])
		categorization[actualAve == 0] = -2
		categorization[actualAve == targetAve] = -1
		for i, threshold in enumerate(thresholds):
			# categorization[targetAve / actualAve > threshold] = i + 1
			categorization[actualAve / targetAve > threshold] = i + 1

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
		output.writerow(["Km and kcat", "Target", "Actual", "Category"])
		for reaction, target, flux, category in zip(constrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], categorization[kmAndKcatReactions]):
			output.writerow([reaction, target, flux, category])

		output.writerow(["kcat only"])
		for reaction, target, flux, category in zip(constrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], categorization[kcatOnlyReactions]):
			output.writerow([reaction, target, flux, category])

		if np.sum(mixedReactions):
			output.writerow(["mixed constraints"])
			for reaction, target, flux, category in zip(constrainedReactions[mixedReactions], targetAve[mixedReactions], actualAve[mixedReactions], categorization[mixedReactions]):
				output.writerow([reaction, target, flux, category])

		csvFile.close()

		targetAve += 1e-6
		actualAve += 1e-6

		axes_limits = [1e-7, 1e4]
		plt.figure(figsize = (8, 8))
		ax = plt.axes()
		plt.loglog(axes_limits, axes_limits, 'k')
		plt.loglog(targetAve, actualAve, "ob", markeredgewidth = 0.25, alpha = 0.25)
		plt.xlabel("Target Flux (mmol/g/hr)")
		plt.ylabel("Actual Flux (mmol/g/hr)")
		plt.minorticks_off()
		whitePadSparklineAxis(ax)
		ax.set_ylim(axes_limits)
		ax.set_xlim(axes_limits)
		ax.set_yticks(axes_limits)
		ax.set_xticks(axes_limits)

		exportFigure(plt, plotOutDir, plotOutFileName)
		plt.close("all")

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
			"reset",
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
			"reset",
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

		p = bokeh.io.vplot(p1, p2)
		bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + ".html"), title=plotOutFileName, autosave=False)
		bokeh.io.save(p)
		bokeh.io.curstate().reset()

if __name__ == "__main__":
	Plot().cli()
