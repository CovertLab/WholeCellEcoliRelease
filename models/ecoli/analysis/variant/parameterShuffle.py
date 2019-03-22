from __future__ import absolute_import

import cPickle

import os

from xlrd import open_workbook
import numpy as np
import matplotlib.pyplot as plt

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

TITLE_FONTSIZE = 10
LABEL_FONTSIZE = 8.5
MARKERSIZE = 30
SCATTER_MSIZE = 5
LINEWIDTH = 4

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

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
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
		flag1 = 0  # failure
		flag2 = 0  # 3-hour upper limit
		flag3 = 0  # final dry mass < 750
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
		pfAx = plt.subplot2grid((nrows, ncols), (4, 3 + 1), rowspan = 2, colspan = 2)
		fdAx = plt.subplot2grid((nrows, ncols), (4, 6 + 1), rowspan = 2, colspan = 2)
		dp_dAx = plt.subplot2grid((nrows, ncols), (3, 1), rowspan = 1, colspan = 2)
		pf_pAx = plt.subplot2grid((nrows, ncols), (3, 4), rowspan = 1, colspan = 2)
		dp_pAx = plt.subplot2grid((nrows, ncols), (4, 0), rowspan = 2, colspan = 1)
		fd_fAx = plt.subplot2grid((nrows, ncols), (3, 7), rowspan = 1, colspan = 2)
		pf_fAx = plt.subplot2grid((nrows, ncols), (4, 3), rowspan = 2, colspan = 1)
		fd_dAx = plt.subplot2grid((nrows, ncols), (4, 6), rowspan = 2, colspan = 1)
		histAxes = [diviAx, protAx, fluxAx]
		axesList = [diviAx, protAx, fluxAx, dpAx, pfAx, fdAx, dp_dAx, dp_pAx, pf_pAx, fd_fAx, pf_fAx, fd_dAx]

		# Plot histograms
		nBins = 50
		histTitles = ["Division Time", "Proteome Comparison", "Fluxome Comparison"]
		histXlabels = ["Time (min)", "PCC", "PCC"]
		for i, ax in enumerate(histAxes):
			if i == 2:
				fluxAx.set_ylim([fluxAx.get_ylim()[0], 55])
			ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
			ax.plot([wt[i], wt[i]], [0, ax.get_ylim()[1]], "r", linewidth = LINEWIDTH)
			ax.set_title("%s\nwt = %0.2f" % (histTitles[i], wt[i]), fontsize = TITLE_FONTSIZE)
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
		axesList = [dpAx, pfAx, fdAx]
		xIndices = [0, 1, 2]
		yIndices = [1, 2, 0]
		for i, ax in enumerate(axesList):
			ax.scatter(data[:, xIndices[i]], data[:, yIndices[i]], s = 2, edgecolor = "none", facecolor = "k")
			ax.plot(wt[xIndices[i]], wt[yIndices[i]], "r.", markersize = MARKERSIZE)
			ax.set_xlim(histAxes[xIndices[i]].get_xlim())
			ax.set_ylim(histAxes[yIndices[i]].get_xlim())
			whitePadSparklineAxis(ax)
			ax.axvline(x = wt[xIndices[i]], linestyle = "--", color = "r")
			ax.axhline(y = wt[yIndices[i]], linestyle = "--", color = "r")

		# Plot small histograms
		fd_fAx.set_ylim([fluxAx.get_ylim()[0], 55])
		titles = ["Division Time (min)", "Proteome Comparison (PCC)", "Fluxome Comparison (PCC)"]
		for i, ax in enumerate([dp_dAx, pf_pAx, fd_fAx]):
			ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
			ax.axvline(wt[i], color = "r", linewidth = LINEWIDTH)
			removeXAxis(ax)
			whitePadSparklineAxis(ax)
			ax.set_title(titles[i], fontsize = TITLE_FONTSIZE)

		pf_fAx.set_xlim([fluxAx.get_ylim()[0], 55])
		for i, ax in enumerate([fd_dAx, dp_pAx, pf_fAx]):
			ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5, orientation = "horizontal")
			ax.axhline(wt[i], color = "r", linewidth = LINEWIDTH)
			removeYAxis(ax)
			whitePadSparklineAxis(ax)
			ax.set_ylabel(titles[i], fontsize = TITLE_FONTSIZE)
			ax.invert_xaxis()

		for ax in axesList:
			ax.tick_params(labelsize = LABEL_FONTSIZE)

		plt.subplots_adjust(wspace = 2, hspace = 2, left = 0.075, right = 0.95)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		if GENERATE_CLEAN_SUBPLOTS:

			# Histograms
			fig, ax = plt.subplots(1, 1, figsize = (5, 5))
			filenames = ["Division", "Proteome", "Fluxome"]
			axesList = [diviAx, protAx, fluxAx]
			for i, filename in enumerate(filenames):
				ax.set_xlim(axesList[i].get_xlim())
				ax.set_ylim(axesList[i].get_ylim())
				ax.plot([wt[i], wt[i]], [0, ax.get_ylim()[1]], "r", linewidth = LINEWIDTH)
				ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5)
				whitePadSparklineAxis(ax)
				removeXAxis(ax)
				cleanAxis(ax)
				exportFigure(plt, plotOutDir, plotOutFileName + "_%s" % filename, metadata)
				plt.cla()

			# Scatter plots
			fig, ax = plt.subplots(1, 1, figsize = (5, 5))
			filenames = ["DP", "PF", "FD"]
			axesList = [dpAx, pfAx, fdAx]
			xIndices = [0, 1, 2]
			yIndices = [1, 2, 0]
			for i, filename in enumerate(filenames):
				ax.scatter(data[:, xIndices[i]], data[:, yIndices[i]], s = SCATTER_MSIZE, edgecolor = "none", facecolor = "k")
				ax.plot(wt[xIndices[i]], wt[yIndices[i]], "r.", markersize = MARKERSIZE)
				ax.set_xlim(axesList[i].get_xlim())
				ax.set_ylim(axesList[i].get_ylim())
				removeYAxis(ax)
				whitePadSparklineAxis(ax)
				ax.axvline(x = wt[xIndices[i]], linestyle = "--", color = "r")
				ax.axhline(y = wt[yIndices[i]], linestyle = "--", color = "r")
				cleanAxis(ax)
				exportFigure(plt, plotOutDir, plotOutFileName + "_%s" % filename, metadata)
				plt.cla()
			plt.close("all")

			# Small histograms for splom (only generates side histograms)
			filenames = ["Division", "Proteome", "Fluxome"]
			axesList = [diviAx, protAx, fluxAx]
			for i, filename in enumerate(filenames):
				fig, ax = plt.subplots(1, 1, figsize = (2, 5))
				ax.set_ylim(axesList[i].get_xlim())
				ax.set_xlim(axesList[i].get_ylim())
				ax.hist(data[:, i], bins = nBins, edgecolor = "none", facecolor = "k", alpha = 0.5, orientation = "horizontal")
				ax.axhline(wt[i], color = "r", linewidth = LINEWIDTH)
				ax.invert_xaxis()
				whitePadSparklineAxis(ax)
				cleanAxis(ax)
				exportFigure(plt, plotOutDir, plotOutFileName + "_%s_%s" % (filename, "small"), metadata)
				plt.close("all")

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
