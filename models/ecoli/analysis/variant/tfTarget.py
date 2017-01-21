import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle
import scipy.stats

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):
	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"
	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)
	variants = sorted(ap._path_data['variant'].tolist()) # Sorry for accessing private data
	variant = variants[0]
	sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))

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
	
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

	# Bokeh
	from bokeh.plotting import figure, output_file, ColumnDataSource, save
	from bokeh.models import HoverTool, BoxZoomTool, LassoSelectTool, PanTool, WheelZoomTool, ResizeTool, UndoTool, RedoTool

	source = ColumnDataSource(data = dict(x = x, y = y[sortedIdxs], targetId = targetIds, tfId = tfs, condition = conditions))
	hover = HoverTool(tooltips = [("target", "@targetId"), ("TF", "@tfId"), ("condition", "@condition")])
	tools = [hover, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]
	plot = figure(x_axis_label = xlabel, y_axis_label = ylabel, width = 800, height = 500, tools = tools)

	plot.scatter("x", "y", source = source)

	if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
		os.makedirs(os.path.join(plotOutDir, "html_plots"))
	import bokeh.io
	bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + "__probBound" + ".html"), title = plotOutFileName, autosave = False)
	bokeh.io.save(plot)

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