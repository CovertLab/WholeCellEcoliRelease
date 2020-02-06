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

import numpy as np
from matplotlib import pyplot as plt
import bokeh.io
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import (HoverTool, BoxZoomTool, LassoSelectTool, PanTool,
	WheelZoomTool, ResizeTool, UndoTool, RedoTool)

from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

BURN_IN_STEPS = 20


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# read constraint data
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		allTargetFluxes = enzymeKineticsReader.readColumn("targetFluxes")
		allActualFluxes = enzymeKineticsReader.readColumn("actualFluxes")
		kineticsConstrainedReactions = np.array(enzymeKineticsReader.readAttribute("kineticsConstrainedReactions"))
		constraint_is_kcat_only = np.array(enzymeKineticsReader.readAttribute('constraint_is_kcat_only'))

		# kinetic target fluxes
		targetFluxes = allTargetFluxes[:, 0:len(kineticsConstrainedReactions)]
		actualFluxes = allActualFluxes[:, 0:len(kineticsConstrainedReactions)]

		targetAve = np.mean(targetFluxes[BURN_IN_STEPS:, :], axis = 0)
		actualAve = np.mean(actualFluxes[BURN_IN_STEPS:, :], axis = 0)

		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only

		kmAndKcatThresholds = [2, 10]
		kmAndKcatCategorization = np.zeros(np.sum(kmAndKcatReactions))
		for i, threshold in enumerate(kmAndKcatThresholds):
			kmAndKcatCategorization[targetAve[kmAndKcatReactions] / actualAve[kmAndKcatReactions] > threshold] = i + 1
			kmAndKcatCategorization[actualAve[kmAndKcatReactions] / targetAve[kmAndKcatReactions] > threshold] = i + 1
		kmAndKcatCategorization[actualAve[kmAndKcatReactions] == 0] = -1

		kcatOnlyThresholds = [2, 10]
		kcatOnlyCategorization = np.zeros(np.sum(kcatOnlyReactions))
		for i, threshold in enumerate(kcatOnlyThresholds):
			kcatOnlyCategorization[actualAve[kcatOnlyReactions] / targetAve[kcatOnlyReactions] > threshold] = i + 1
		kcatOnlyCategorization[actualAve[kcatOnlyReactions] == 0] = -1

		targetAve += 1e-13
		actualAve += 1e-13

		plt.figure(figsize = (8, 8))
		plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 0], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 0], "og")
		plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 1], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 1], "o")
		plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 2], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 2], "or")
		plt.loglog([1e-12, 1], [1e-12, 1], '--g')
		plt.loglog([1e-12, 1], [1e-11, 10], '--r')
		plt.xlabel("Target Flux (dmol/L/s)")
		plt.ylabel("Actual Flux (dmol/L/s)")

		exportFigure(plt, plotOutDir, plotOutFileName)
		plt.close("all")

		return

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
