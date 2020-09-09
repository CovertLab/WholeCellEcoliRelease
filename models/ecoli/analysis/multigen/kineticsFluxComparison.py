"""
Compare fluxes in simulation to target fluxes

@date: Created 12/15/16
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import io
import os
import re

import bokeh.io
import bokeh.io.state
import bokeh.layouts
from bokeh.models import HoverTool
from bokeh.plotting import figure, ColumnDataSource
from matplotlib import pyplot as plt
import numpy as np
from six.moves import cPickle, zip

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.io.tablereader import TableReader
from wholecell.io import tsv
from wholecell.utils import filepath, units
from wholecell.utils.sparkline import whitePadSparklineAxis


BURN_IN_STEPS = 20


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()
		# allDir = ap.get_cells(generation = [0, 1, 2])

		sim_data = cPickle.load(open(simDataFile, "rb"))

		allTargetFluxList = []
		allActualFluxList = []

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			mainListener = TableReader(os.path.join(simOutDir, "Main"))
			initialTime = mainListener.readAttribute("initialTime")
			time = mainListener.readColumn("time") - initialTime
			mainListener.close()

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			massListener.close()

			coefficient = dryMass / cellMass * sim_data.constants.cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

			# read constraint data
			enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
			allTargetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
			allActualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
			kineticsConstrainedReactions = np.array(enzymeKineticsReader.readAttribute("kineticsConstrainedReactions"))
			enzymeKineticsReader.close()

			allTargetFluxes = allTargetFluxes.asNumber(units.mmol / units.g / units.h)
			allActualFluxes = allActualFluxes.asNumber(units.mmol / units.g / units.h)

			allTargetAve = np.mean(allTargetFluxes[BURN_IN_STEPS:, :], axis=0)
			allActualAve = np.mean(allActualFluxes[BURN_IN_STEPS:, :], axis=0)

			if len(allTargetFluxList) == 0:
				allTargetFluxList = np.array([allTargetAve])
				allActualFluxList = np.array([allActualAve])
			else:
				allTargetFluxList = np.concatenate((allTargetFluxList, np.array([allTargetAve])), axis = 0)
				allActualFluxList = np.concatenate((allActualFluxList, np.array([allActualAve])), axis = 0)

		allTargetAve = np.mean(allTargetFluxList, axis = 0)
		allActualAve = np.mean(allActualFluxList, axis = 0)

		n_kinetic_constrained_reactions = len(kineticsConstrainedReactions)

		# boundary target fluxes
		boundaryTargetAve = allTargetAve[n_kinetic_constrained_reactions:]
		boundaryActualAve = allActualAve[n_kinetic_constrained_reactions:]

		# kinetic target fluxes
		actualFluxes = allActualFluxes[:, :n_kinetic_constrained_reactions]
		targetAve = allTargetAve[:n_kinetic_constrained_reactions]
		actualAve = allActualAve[:n_kinetic_constrained_reactions]

		thresholds = [2, 10]
		categorization = np.zeros(n_kinetic_constrained_reactions)
		categorization[actualAve == 0] = -2
		categorization[actualAve == targetAve] = -1
		for i, threshold in enumerate(thresholds):
			# categorization[targetAve / actualAve > threshold] = i + 1
			categorization[actualAve / targetAve > threshold] = i + 1

		# url for ecocyc to highlight fluxes that are 0 on metabolic network diagram
		siteStr = "https://ecocyc.org/overviewsWeb/celOv.shtml?zoomlevel=1&orgid=ECOLI"
		excluded = ['RXN0-2201', 'RXN-16000', 'RXN-12583', 'RXN-11496', 'DIMESULFREDUCT-RXN', '3.6.1.41-R[4/63051]5-NUCLEOTID-RXN'] # reactions not recognized by ecocyc
		rxns = []
		for i, reaction in enumerate(kineticsConstrainedReactions[categorization == -2]):
			rxn = re.findall(".+RXN", reaction)
			if len(rxn) == 0:
				rxn = re.findall("RXN[^-]*-[0-9]+", reaction)
			if rxn[0] not in excluded:
				siteStr += "&rnids=%s" % rxn[0]
			rxns.append(rxn[0])
		# print(siteStr)

		csvFile = io.open(os.path.join(plotOutDir, plotOutFileName + ".tsv"), "wb")
		output = tsv.writer(csvFile)
		output.writerow(["ecocyc link:", siteStr])
		output.writerow(["", "Target", "Actual", "Category"])

		for reaction, target, flux, category in zip(kineticsConstrainedReactions, targetAve, actualAve, categorization):
			output.writerow([reaction, target, flux, category])

		csvFile.close()

		targetAve += 1e-6
		actualAve += 1e-6

		axes_limits = [1e-7, 1e4]
		plt.figure(figsize = (8, 8))
		ax = plt.axes()
		plt.loglog(axes_limits, axes_limits, 'k')
		plt.loglog(targetAve, actualAve, "ob", markeredgewidth = 0.25, alpha = 0.25)
		plt.loglog(boundaryTargetAve, boundaryActualAve, "ob", c='r', markeredgewidth=0.25, alpha=0.9, label='boundary fluxes')
		plt.xlabel("Target Flux (mmol/g/hr)")
		plt.ylabel("Actual Flux (mmol/g/hr)")
		plt.minorticks_off()
		whitePadSparklineAxis(ax)
		# noinspection PyTypeChecker
		ax.set_ylim(axes_limits)
		# noinspection PyTypeChecker
		ax.set_xlim(axes_limits)
		ax.set_yticks(axes_limits)
		ax.set_xticks(axes_limits)
		ax.legend()

		exportFigure(plt, plotOutDir, plotOutFileName)
		plt.close("all")

		source = ColumnDataSource(
			data = dict(
				x = targetAve,
				y = actualAve,
				reactionName = kineticsConstrainedReactions)
			)
		hover = HoverTool(
			tooltips = [
				("Reaction", "@reactionName"),
				]
			)
		tools = [
			hover,
			'box_zoom',
			'lasso_select',
			'pan',
			'wheel_zoom',
			'undo',
			'redo',
			'reset',
			]
		p1 = figure(
			x_axis_label="Target",
			x_axis_type="log",
			x_range=[min(targetAve[targetAve > 0]), max(targetAve)],
			y_axis_label="Actual",
			y_axis_type="log",
			y_range=[min(actualAve[actualAve > 0]), max(actualAve)],
			width=800,
			height=800,
			tools=tools,
			)
		p1.scatter('x', 'y', source=source, size=8)
		p1.line([1e-15, 10], [1e-15, 10], line_color="red", line_dash="dashed")

		## flux for each reaction
		hover2 = HoverTool(
			tooltips = [
				("Reaction", "@reactionName"),
				]
			)
		tools2 = [
			hover2,
			'box_zoom',
			'lasso_select',
			'pan',
			'wheel_zoom',
			'undo',
			'redo',
			'reset',
			]
		p2 = figure(
			x_axis_label="Time(s)",
			y_axis_label="Flux",
			y_axis_type="log",
			y_range=[1e-8, 1],
			width=800,
			height=800,
			tools=tools2,
			)

		colors = COLORS_LARGE
		nTimesteps = len(time[BURN_IN_STEPS:])
		x = time[BURN_IN_STEPS:]
		y = actualFluxes[BURN_IN_STEPS:, 0]
		reactionName = np.repeat(kineticsConstrainedReactions[0], nTimesteps)

		source = ColumnDataSource(
			data = dict(
				x = x,
				y = y,
				reactionName = reactionName)
			)
		p2.line('x', 'y', line_color=colors[0], source=source)

		# Plot remaining metabolites onto initialized figure
		for m in np.arange(1, actualFluxes.shape[1]):
			y = actualFluxes[BURN_IN_STEPS:, m]
			reactionName = np.repeat(kineticsConstrainedReactions[m], nTimesteps)

			source = ColumnDataSource(
				data = dict(
					x = x,
					y = y,
					reactionName = reactionName)
			)
			p2.line('x', 'y', line_color=colors[m % len(colors)], source=source)

		html_dir = filepath.makedirs(plotOutDir, "html_plots")
		p = bokeh.layouts.gridplot([[p1], [p2]])
		bokeh.io.output_file(os.path.join(html_dir, plotOutFileName + ".html"), title=plotOutFileName)
		bokeh.io.save(p)
		bokeh.io.state.curstate().reset()


if __name__ == "__main__":
	Plot().cli()
