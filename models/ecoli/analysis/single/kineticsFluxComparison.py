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

from models.ecoli.analysis import singleAnalysisPlot
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.io.tablereader import TableReader
from wholecell.io import tsv
from wholecell.utils import filepath, units
from wholecell.utils.sparkline import whitePadSparklineAxis


BURN_IN_STEPS = 20


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, 'rb'))

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
		constraint_is_kcat_only = np.array(enzymeKineticsReader.readAttribute('constraint_is_kcat_only'))

		allTargetFluxes = allTargetFluxes.asNumber(units.mmol / units.g / units.h)
		allActualFluxes = allActualFluxes.asNumber(units.mmol / units.g / units.h)

		allTargetAve = np.mean(allTargetFluxes[BURN_IN_STEPS:, :], axis = 0)
		allActualAve = np.mean(allActualFluxes[BURN_IN_STEPS:, :], axis = 0)

		n_kinetic_constrained_reactions = len(kineticsConstrainedReactions)

		# boundary target fluxes
		boundaryTargetAve = allTargetAve[n_kinetic_constrained_reactions:]
		boundaryActualAve = allActualAve[n_kinetic_constrained_reactions:]

		# kinetic target fluxes
		actualFluxes = allActualFluxes[:, :n_kinetic_constrained_reactions]
		targetAve = allTargetAve[:n_kinetic_constrained_reactions]
		actualAve = allActualAve[:n_kinetic_constrained_reactions]

		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only

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
		for i, reaction in enumerate(kineticsConstrainedReactions):
			if actualAve[i] == 0:
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
		output.writerow(["Km and kcat", "Target", "Actual", "Category"])
		for reaction, target, flux, category in zip(kineticsConstrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], categorization[kmAndKcatReactions]):
			output.writerow([reaction, target, flux, category])

		output.writerow(["kcat only"])
		for reaction, target, flux, category in zip(kineticsConstrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], categorization[kcatOnlyReactions]):
			output.writerow([reaction, target, flux, category])

		csvFile.close()

		# Masks for highlighting and counting certain groups
		on_target_mask = actualAve == targetAve
		zero_mask = (actualAve == 0) & ~on_target_mask
		other_mask = ~on_target_mask & ~zero_mask

		targetAve += 1e-6
		actualAve += 1e-6

		axes_limits = [1e-7, 1e4]
		plt.figure(figsize = (8, 8))
		ax = plt.axes()
		plt.loglog(axes_limits, axes_limits, 'k--')
		plt.loglog(targetAve[on_target_mask], actualAve[on_target_mask], 'og',
			markeredgewidth=0.25, alpha=0.25,
			label='on target fluxes (n={})'.format(on_target_mask.sum()))
		plt.loglog(targetAve[zero_mask], actualAve[zero_mask], 'or',
			markeredgewidth=0.25, alpha=0.25,
			label='zero fluxes (n={})'.format(zero_mask.sum()))
		plt.loglog(targetAve[other_mask], actualAve[other_mask], 'ob',
			markeredgewidth=0.25, alpha=0.25,
			label='other fluxes (n={})'.format(other_mask.sum()))
		if len(boundaryTargetAve):
			plt.loglog(boundaryTargetAve, boundaryActualAve, 'ok',
				markeredgewidth=0.25, alpha=0.9, label='boundary fluxes')
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
