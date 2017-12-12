from xlrd import open_workbook
import numpy as np
import matplotlib.pyplot as plt
import cPickle
from wholecell.utils.sparkline import whitePadSparklineAxis
import wholecell.utils.constants
from wholecell.utils import units
import argparse
import os
import re

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

TITLE_FONTSIZE = 10
LABEL_FONTSIZE = 8.5
MARKERSIZE = 20
SCATTER_MSIZE = 5

GENERATE_CLEAN_SUBPLOTS = True

def removeXAxis(ax):
	ax.spines["bottom"].set_visible(False)
	ax.tick_params(bottom = "off")
	ax.set_xticklabels([])

def removeYAxis(ax):
	ax.spines["left"].set_visible(False)
	ax.tick_params(left = "off")
	ax.set_yticklabels([])

def cleanAxis(ax):
	ax.set_xticklabels([])
	ax.set_yticklabels([])

def returnAxis(ax):
	ax.spines["bottom"].set_visible(True)
	ax.spines["left"].set_visible(True)
	ax.tick_params(bottom = "on", left = "on")

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):
	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	if not os.path.exists(os.path.join(plotOutDir, "distribution_division_fluxome_proteome_data_matrix.xls")):
		print "%s must exist as an .xls file for this plot." % "distribution_division_fluxome_proteome_data_matrix.xls"
		return	

	# Load data
	wb = open_workbook(os.path.join(plotOutDir, "distribution_division_fluxome_proteome_data_matrix.xls"))
	data_raw = wb.sheet_by_index(0)
	heading = data_raw.row_values(0)
	col_divisionTime = heading.index("division time")
	col_proteome = heading.index("proteome pearson r")
	col_fluxome = heading.index("fluxome pearson r")
	col_finalMass = heading.index("final dry mass (fg)")
	dataIndices = [col_divisionTime, col_proteome, col_fluxome]
	wt = np.array(data_raw.row_values(1), dtype = float)[dataIndices]

	# Load failures
	failures = cPickle.load(open(os.path.join(plotOutDir, "failed_variants.cPickle"), "rb"))
	flag1 = 0; # failure
	flag2 = 0; # 3-hour upper limit
	flag3 = 0; # final dry mass < 750
	data = []
	for varId, row in enumerate(xrange(2, data_raw.nrows)):
		row_values = data_raw.row_values(row)
		
		if (varId + 1) in failures:
			flag1 += 1
		elif row_values[col_divisionTime] == 180:
			flag2 += 1
		elif row_values[col_finalMass] < 750.:
			flag3 += 1
		else:
			data.append(np.array(row_values, dtype = float)[dataIndices])
	data = np.array(data)

	# Plot
	nrows = 6
	ncols = 9
	factor = 1.5
	fig = plt.figure(figsize = (ncols * factor, nrows * factor))
	diviAx = plt.subplot2grid((nrows, ncols), (0, 0 + 1), rowspan = 2, colspan = 2)
	protAx = plt.subplot2grid((nrows, ncols), (0, 3 + 1), rowspan = 2, colspan = 2)
	fluxAx = plt.subplot2grid((nrows, ncols), (0, 6 + 1), rowspan = 2, colspan = 2)
	dpAx = plt.subplot2grid((nrows, ncols), (4, 0 + 1), rowspan = 2, colspan = 2)
	dfAx = plt.subplot2grid((nrows, ncols), (4, 3 + 1), rowspan = 2, colspan = 2)
	pfAx = plt.subplot2grid((nrows, ncols), (4, 6 + 1), rowspan = 2, colspan = 2)
	dp_dAx = plt.subplot2grid((nrows, ncols), (3, 1), rowspan = 1, colspan = 2)
	df_dAx = plt.subplot2grid((nrows, ncols), (3, 4), rowspan = 1, colspan = 2)
	dp_pAx = plt.subplot2grid((nrows, ncols), (4, 0), rowspan = 2, colspan = 1)
	pf_pAx = plt.subplot2grid((nrows, ncols), (3, 7), rowspan = 1, colspan = 2)
	df_fAx = plt.subplot2grid((nrows, ncols), (4, 3), rowspan = 2, colspan = 1)
	pf_fAx = plt.subplot2grid((nrows, ncols), (4, 6), rowspan = 2, colspan = 1)
	histAxes = [diviAx, protAx, fluxAx]
	axesList = [diviAx, protAx, fluxAx, dpAx, dfAx, pfAx, dp_dAx, dp_pAx, df_dAx, pf_pAx, df_fAx, pf_fAx]

	# Plot histograms
	nBins = 50
	histTitles = ["Division Time", "Proteome Comparison", "Fluxome Comparison"]
	histXlabels = ["Time (min)", "PCC", "PCC"]
	for i, ax in enumerate(histAxes):
		if i == 2:
			fluxAx.set_ylim([fluxAx.get_ylim()[0], 55])
		ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
		ax.plot([wt[i], wt[i]], [0, ax.get_ylim()[1]], "b")
		ax.set_title(histTitles[i], fontsize = TITLE_FONTSIZE)
		ax.set_ylabel("Frequency", fontsize = TITLE_FONTSIZE)
		ax.set_xlabel(histXlabels[i], fontsize = TITLE_FONTSIZE)

	diviAx.set_xlim([30, diviAx.get_xlim()[1]])
	diviAx.text(wt[0] - 5, 50, "%0.0f %%" % ((data[:, 0] < wt[0]).sum() / data.shape[0] * 100.), size = LABEL_FONTSIZE, horizontalalignment = "right")
	diviAx.text(wt[0] + 5, 50, "%0.0f %%" % ((data[:, 0] > wt[0]).sum() / data.shape[0] * 100.), size = LABEL_FONTSIZE, horizontalalignment = "left")

	protAx.text(wt[1] - 0.005, 40, "%0.0f %%" % ((data[:, 1] < wt[1]).sum() / data.shape[0] * 100.), size = LABEL_FONTSIZE, horizontalalignment = "right")
	protAx.text(wt[1] + 0.005, 40, "%0.0f %%" % ((data[:, 1] > wt[1]).sum() / data.shape[0] * 100.), size = LABEL_FONTSIZE, horizontalalignment = "left")

	fluxAx.text(wt[2] - 0.1, 50, "%0.1f %%" % (float((data[:, 2] < wt[2]).sum()) / data.shape[0] * 100.), size = LABEL_FONTSIZE, horizontalalignment = "right")
	fluxAx.text(wt[2] + 0.1, 50, "%0.1f %%" % (float((data[:, 2] > wt[2]).sum()) / data.shape[0] * 100.), size = LABEL_FONTSIZE, horizontalalignment = "left")
	
	for ax in histAxes:
		whitePadSparklineAxis(ax)

	# Plot scatter plots
	dpAx.scatter(data[:, 0], data[:, 1], s = 2, edgecolor = "none", facecolor = "k")
	dpAx.plot(wt[0], wt[1], "b.", markersize = MARKERSIZE)
	dpAx.set_xlim(diviAx.get_xlim())
	dpAx.set_ylim(protAx.get_xlim())
	whitePadSparklineAxis(dpAx)
	dpAx.axvline(x = wt[0], linestyle = "--")
	dpAx.axhline(y = wt[1], linestyle = "--")

	dfAx.scatter(data[:, 0], data[:, 2], s = 2, edgecolor = "none", facecolor = "k")
	dfAx.plot(wt[0], wt[2], "b.", markersize = MARKERSIZE)
	dfAx.set_xlim(diviAx.get_xlim())
	dfAx.set_ylim(fluxAx.get_xlim())
	whitePadSparklineAxis(dfAx)
	dfAx.axvline(x = wt[0], linestyle = "--")
	dfAx.axhline(y = wt[2], linestyle = "--")

	pfAx.scatter(data[:, 1], data[:, 2], s = 2, edgecolor = "none", facecolor = "k")
	pfAx.plot(wt[1], wt[2], "b.", markersize = MARKERSIZE)
	pfAx.set_xlim(protAx.get_xlim())
	pfAx.set_ylim(fluxAx.get_xlim())
	whitePadSparklineAxis(pfAx)
	pfAx.axvline(x = wt[1], linestyle = "--")
	pfAx.axhline(y = wt[2], linestyle = "--")

	# Plot small histograms
	for ax in [dp_dAx, df_dAx]:
		ax.hist(data[:, 0], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
		ax.plot([wt[0], wt[0]], [0, ax.get_ylim()[1]], "b")
		removeXAxis(ax)
		whitePadSparklineAxis(ax)

	for ax in [dp_pAx, pf_pAx]:
		orientation = "vertical"
		if ax == dp_pAx:
			orientation = "horizontal"
		ax.hist(data[:, 1], bins = nBins, edgecolor = "none", orientation = orientation, facecolor = "k", alpha = 0.5)
		if orientation == "horizontal":
			ax.plot([0, ax.get_xlim()[1]], [wt[1], wt[1]], "b")
			removeYAxis(ax)
		else:
			ax.plot([wt[1], wt[1]], [0, ax.get_ylim()[1]], "b")
			removeXAxis(ax)
		whitePadSparklineAxis(ax)

	for ax in [df_fAx, pf_fAx]:
		ax.hist(data[:, 2], bins = nBins, edgecolor = "none", orientation = "horizontal", facecolor = "k", alpha = 0.5)
		ax.plot([0, ax.get_xlim()[1]], [wt[2], wt[2]], "b")
		ax.set_xlim([0, 55])
		removeYAxis(ax)
		whitePadSparklineAxis(ax)

	dp_dAx.set_title("Division Time (min)", fontsize = TITLE_FONTSIZE)
	df_dAx.set_title("Division Time (min)", fontsize = TITLE_FONTSIZE)
	dp_pAx.set_ylabel("Proteome Comparison (PCC)", fontsize = TITLE_FONTSIZE, rotation = "vertical")
	pf_pAx.set_title("Proteome Comparison (PCC)", fontsize = TITLE_FONTSIZE)
	df_fAx.set_ylabel("Fluxome Comparison (PCC)", fontsize = TITLE_FONTSIZE, rotation = "vertical")
	pf_fAx.set_ylabel("Fluxome Comparison (PCC)", fontsize = TITLE_FONTSIZE, rotation = "vertical")

	for ax in [dp_pAx, df_fAx, pf_fAx]:
		ax.invert_xaxis()
		ax.xaxis.set_label_position("top")

	for ax in axesList:
		ax.tick_params(labelsize = LABEL_FONTSIZE)

	plt.subplots_adjust(wspace = 2, hspace = 2, left = 0.075, right = 0.95)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	if GENERATE_CLEAN_SUBPLOTS:

		# Histograms
		fig, ax = plt.subplots(1, 1, figsize = (6, 6))
		filenames = ["Division", "Proteome", "Fluxome"]
		axesList = [diviAx, protAx, fluxAx]
		for i, filename in enumerate(filenames):
			ax.set_xlim(axesList[i].get_xlim())
			ax.set_ylim(axesList[i].get_ylim())
			ax.plot([wt[i], wt[i]], [0, ax.get_ylim()[1]], "b")
			ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
			whitePadSparklineAxis(ax)
			cleanAxis(ax)
			exportFigure(plt, plotOutDir, plotOutFileName + "_%s" % filename, metadata)
			plt.cla()

		# Sploms
		fig, ax = plt.subplots(1, 1, figsize = (5, 5))
		filenames = ["DP", "DF", "PF"]
		axesList = [dpAx, dfAx, pfAx]
		xIndices = [0, 0, 1]
		yIndices = [1, 2, 2]
		for i, filename in enumerate(filenames):
			ax.scatter(data[:, xIndices[i]], data[:, yIndices[i]], s = SCATTER_MSIZE, edgecolor = "none", facecolor = "k")
			ax.plot(wt[xIndices[i]], wt[yIndices[i]], "b.", markersize = MARKERSIZE)
			ax.set_xlim(axesList[i].get_xlim())
			ax.set_ylim(axesList[i].get_ylim())
			removeYAxis(ax)
			whitePadSparklineAxis(ax)
			ax.axvline(x = wt[xIndices[i]], linestyle = "--")
			ax.axhline(y = wt[yIndices[i]], linestyle = "--")
			cleanAxis(ax)
			exportFigure(plt, plotOutDir, plotOutFileName + "_%s" % filename, metadata)
			plt.cla()
		plt.close("all")

		# Small histograms for splom
		filenames = ["Division", "Proteome", "Fluxome"]
		axesList = [diviAx, protAx, fluxAx]
		for i, filename in enumerate(filenames):

			if i in [0, 1]:	
				fig, ax = plt.subplots(1, 1, figsize = (5, 2))
				ax.set_xlim(axesList[i].get_xlim())
				ax.set_ylim(axesList[i].get_ylim())
				ax.plot([wt[i], wt[i]], [0, ax.get_ylim()[1]], "b")
				ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
				removeXAxis(ax)
				whitePadSparklineAxis(ax)
				cleanAxis(ax)
				exportFigure(plt, plotOutDir, plotOutFileName + "_%s_%s_vertical" % (filename, "small"), metadata)
				plt.close("all")

			if i in [1, 2]:
				fig, ax = plt.subplots(1, 1, figsize = (2, 5))
				ax.set_ylim(axesList[i].get_xlim())
				ax.set_xlim(axesList[i].get_ylim())
				ax.plot([0, ax.get_xlim()[1]], [wt[i], wt[i]], "b")				
				ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5, orientation = "horizontal")
				ax.invert_xaxis()
				whitePadSparklineAxis(ax)
				cleanAxis(ax)
				exportFigure(plt, plotOutDir, plotOutFileName + "_%s_%s_horizontal" % (filename, "small"), metadata)
				plt.close("all")

	plt.close("all")

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
