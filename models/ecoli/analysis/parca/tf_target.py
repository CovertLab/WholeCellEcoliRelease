from __future__ import absolute_import, division, print_function

import cPickle
import os

import numpy as np
from matplotlib import pyplot as plt
import bokeh.io
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import (HoverTool, BoxZoomTool, LassoSelectTool, PanTool,
	WheelZoomTool, ResizeTool, UndoTool, RedoTool)

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import filepath


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = cPickle.load(f)

		targetToFC = {}
		targetToFCTF = {}

		for tf in sim_data.tfToActiveInactiveConds:
			for target in sim_data.tfToFC[tf]:
				if target not in targetToFC:
					targetToFC[target] = []
					targetToFCTF[target] = []
				targetToFC[target].append(np.log2(sim_data.tfToFC[tf][target]))
				targetToFCTF[target].append(tf)

		for target in targetToFC:
			targetToFC[target] = np.array(targetToFC[target])

		targets = sorted(targetToFC)

		x = []
		y = []
		maxVals = []
		tfs = []
		targetIds = []

		for idx, target in enumerate(targets):
			for FC, tf in zip(targetToFC[target], targetToFCTF[target]):
				x.append(idx)
				y.append(FC)

				if targetToFC[target].max() >= -1. * targetToFC[target].min():
					maxVals.append(targetToFC[target].max())
				else:
					maxVals.append(targetToFC[target].min())

				tfs.append(tf)
				targetIds.append(target)
		conditions = [sim_data.conditions[tf + "__active"]["nutrients"] for tf in tfs]

		x = np.array(x)
		y = np.array(y)
		maxVals = np.array(maxVals)

		sortedIdxs = np.argsort(maxVals)
		conditions = [conditions[i] for i in sortedIdxs]
		tfs = [tfs[i] for i in sortedIdxs]
		targetIds = [targetIds[i] for i in sortedIdxs]

		fig = plt.figure(figsize = (11, 8.5))
		ax = plt.subplot(1, 1, 1)
		ax.plot(x, y[sortedIdxs], ".")
		xlabel = "Gene targets (sorted)"
		ylabel = "log2 (Target expression fold change)"
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)

		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close("all")

		source = ColumnDataSource(data = dict(x = x, y = y[sortedIdxs], targetId = targetIds, tfId = tfs, condition = conditions))
		hover = HoverTool(tooltips = [("target", "@targetId"), ("TF", "@tfId"), ("condition", "@condition")])
		tools = [hover, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]
		plot = figure(x_axis_label = xlabel, y_axis_label = ylabel, width = 800, height = 500, tools = tools)

		plot.scatter("x", "y", source = source)

		html_dir = filepath.makedirs(plot_out_dir, 'html_plots')
		html_file = os.path.join(html_dir, plot_out_filename + "__probBound.html")
		bokeh.io.output_file(html_file, title=plot_out_filename, autosave=False)
		bokeh.io.save(plot)
		bokeh.io.curstate().reset()


if __name__ == "__main__":
	Plot().cli()
