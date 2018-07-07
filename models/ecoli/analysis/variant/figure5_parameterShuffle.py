from __future__ import absolute_import

from xlrd import open_workbook
import numpy as np
import matplotlib.pyplot as plt
import cPickle

import os

from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

PLOT_FIG1 = True
PLOT_FIG2 = True
CLEAN_VER = False


def emptyAxis(ax):
	ax.spines["top"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.tick_params(top = "off", bottom = "off", left = "off", right = "off", )
	ax.set_xticks([])
	ax.set_yticks([])

def cleanAxis(ax):
	ax.tick_params(right = "off", top = "off", which = "both", direction = "out", labelsize = 10)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_position(("outward", 10))
	ax.spines["bottom"].set_position(("outward", 10))
	ax.set_yticks(ax.get_ylim())
	ax.set_xticks(ax.get_xlim())

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

		if PLOT_FIG1:
			# Figure 1
			fig1, axesList = plt.subplots(3, 3, figsize = (10, 10))
			allAxesList = axesList.reshape(-1)

			[histDivAx, _, __], [dpAx, histProtAx, ___], [dfAx, pfAx, histFluxAx] = axesList
			histAxes = [histDivAx, histProtAx, histFluxAx]
			histAxesTitles = ["Division Time", "Proteome Comparison", "Fluxome Comparison"]
			histAxesXlabel = ["Time (min)", "PCC", "PCC"]
			emptyAxes = [_, __, ___]

			nBins = 50

			for i, ax in enumerate(histAxes):
				if i == 2:
					histFluxAx.set_ylim([histFluxAx.get_ylim()[0], 55])
				ax.hist(data[:, i], bins = nBins, edgecolor = "none")
				ax.plot([wt[i], wt[i]], [0, ax.get_ylim()[1]], "r")
				if not CLEAN_VER:
					ax.set_title(histAxesTitles[i], y = 1.05)
					ax.set_xlabel(histAxesXlabel[i], fontsize = 10)

			histDivAx.set_xlim([30, histDivAx.get_xlim()[1]])
			if not CLEAN_VER:
				histDivAx.text(wt[0] - 5, 50, "%0.0f %%" % ((data[:, 0] < wt[0]).sum() / data.shape[0] * 100.), size = 8, horizontalalignment = "right")
				histDivAx.text(wt[0] + 5, 50, "%0.0f %%" % ((data[:, 0] > wt[0]).sum() / data.shape[0] * 100.), size = 8, horizontalalignment = "left")

			# histProtAx.set_xlim([histProtAx.get_xlim()[0], 0.76])
			if not CLEAN_VER:
				histProtAx.text(wt[1] - 0.005, 40, "%0.0f %%" % ((data[:, 1] < wt[1]).sum() / data.shape[0] * 100.), size = 8, horizontalalignment = "right")
				histProtAx.text(wt[1] + 0.005, 40, "%0.0f %%" % ((data[:, 1] > wt[1]).sum() / data.shape[0] * 100.), size = 8, horizontalalignment = "left")

			# histFluxAx.set_ylim([histFluxAx.get_ylim()[0], 55])
			if not CLEAN_VER:
				histFluxAx.text(wt[2] - 0.1, 50, "%0.1f %%" % (float((data[:, 2] < wt[2]).sum()) / data.shape[0] * 100.), size = 8, horizontalalignment = "right")
				histFluxAx.text(wt[2] + 0.1, 50, "%0.1f %%" % (float((data[:, 2] > wt[2]).sum()) / data.shape[0] * 100.), size = 8, horizontalalignment = "left")
			# histFluxAx.text(wt[2], 55, "equal: %0.1f %%" % (float((data[:, 2] == wt[2]).sum()) / data.shape[0] * 100.), size = 8, horizontalalignment = "center")

			dpAx.scatter(data[:, 0], data[:, 1], s = 2, edgecolor = "none", facecolor = "k")
			dpAx.plot(wt[0], wt[1], "r.")
			if not CLEAN_VER:
				dpAx.set_xlabel("Division Time (min)", fontsize = 10)
				dpAx.set_ylabel("Proteome Comparison (PCC)", fontsize = 10)
			dpAx.set_xlim(histDivAx.get_xlim())
			dpAx.set_ylim(histProtAx.get_xlim())

			dfAx.scatter(data[:, 0], data[:, 2], s = 2, edgecolor = "none", facecolor = "k")
			dfAx.plot(wt[0], wt[2], "r.")
			if not CLEAN_VER:
				dfAx.set_xlabel("Division Time (min)", fontsize = 10)
				dfAx.set_ylabel("Fluxome Comparison (PCC)", fontsize = 10)
			dfAx.set_xlim(histDivAx.get_xlim())
			dfAx.set_ylim(histFluxAx.get_xlim())

			pfAx.scatter(data[:, 1], data[:, 2], s = 2, edgecolor = "none", facecolor = "k")
			pfAx.plot(wt[1], wt[2], "r.")
			if not CLEAN_VER:
				pfAx.set_xlabel("Proteome Comparison (PCC)", fontsize = 10)
				pfAx.set_ylabel("Fluxome Comparison (PCC)", fontsize = 10)
			pfAx.set_xlim(histProtAx.get_xlim())
			pfAx.set_ylim(histFluxAx.get_xlim())

			for ax in allAxesList:
				cleanAxis(ax)
				if CLEAN_VER:
					ax.set_xticklabels([])
					ax.set_yticklabels([])

			for ax in emptyAxes:
				emptyAxis(ax)

			plt.subplots_adjust(wspace = 0.5, hspace = 0.5)

			if not CLEAN_VER:
				exportFigure(plt, plotOutDir, plotOutFileName + "__splom", metadata)
			else:
				exportFigure(plt, plotOutDir, plotOutFileName + "__splom__clean", metadata)
			plt.close("all")


		if PLOT_FIG2:
		# Figure 2
			fig2 = plt.figure(figsize = (6, 6))
			ax = fig2.add_subplot(1, 1, 1)
			xvals = np.arange(2)
			width = 0.6
			ax.bar(xvals, [1024, flag1 + flag2 + flag3], width, edgecolor = "none")
			if not CLEAN_VER:
				ax.set_xticklabels(["Simulation batch size", "Failures"])
			ax.set_ylim([0, 1024])
			cleanAxis(ax)
			ax.set_xticks(xvals + 0.5 * width)
			ax.set_xlim([ax.get_xticks()[0] - (ax.get_xlim()[1] - ax.get_xticks()[1]), ax.get_xlim()[1]])
			ax.set_yticks([0, flag1 + flag2 + flag3, 1024])

			if not CLEAN_VER:
				exportFigure(plt, plotOutDir, plotOutFileName + "__histogram", metadata)
			else:
				ax.set_yticklabels([])
				ax.set_xticklabels([])
				exportFigure(plt, plotOutDir, plotOutFileName + "__histogram__clean", metadata)
			plt.close("all")


if __name__ == "__main__":
	Plot().cli()
