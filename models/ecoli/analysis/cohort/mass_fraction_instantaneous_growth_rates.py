from __future__ import absolute_import, division, print_function

import os
import numpy as np
from matplotlib import pyplot as plt
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE = 1

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

trim = 0.03
trim_1 = 0.06


def mm2inch(value):
	return value * 0.0393701

def seriesScrubber(series, factor):
	series[abs(series - np.median(series)) > factor * np.nanstd(series)] = np.nan


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(seedOutDir, cohort_plot = True)

		# Get all cells
		allDir = ap.get_cells(seed=[0])

		# fig = plt.figure(figsize = (14, 10))
		# ax1 = plt.subplot2grid((4,2), (0,0), rowspan = 4)
		# ax2 = plt.subplot2grid((4,2), (0,1))
		# ax3 = plt.subplot2grid((4,2), (1,1))
		# ax4 = plt.subplot2grid((4,2), (2,1))
		# ax5 = plt.subplot2grid((4,2), (3,1))
		# axesList = [ax2, ax3, ax4, ax5]
		# colors = ["blue", "green", "red", "cyan"]
		colors = ["#43aa98", "#0071bb", "#bf673c"]

		mult = 2.3

		fig, axis = plt.subplots(1,1)
		fig.set_figwidth(mm2inch(55.2)*mult)
		fig.set_figheight(mm2inch(53.36)*mult)

		width = 200

		timeMultigen = np.zeros(0)
		cellMassGrowthRateMultigen = np.zeros(0)
		proteinGrowthRateMultigen = np.zeros(0)
		rnaGrowthRateMultigen = np.zeros(0)


		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			time = main_reader.readColumn("time")
			timeStep = main_reader.readColumn("timeStepSec")

			mass = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = mass.readColumn("cellMass")
			proteinMass = mass.readColumn("proteinMass")
			rnaMass = mass.readColumn("rnaMass")
			dnaMass = mass.readColumn("dnaMass")

			cellmassGrowthRate = np.diff(cellMass) / cellMass[1:] / timeStep[:-1]
			proteinGrowthRate = np.diff(proteinMass) / proteinMass[1:] / timeStep[:-1]
			rnaGrowthRate = np.diff(rnaMass) / rnaMass[1:] / timeStep[:-1]
			dnaGrowthRate = np.diff(dnaMass) / dnaMass[1:] / timeStep[:-1]

			timeMultigen = np.hstack((timeMultigen, time[:-1]))
			cellMassGrowthRateMultigen = np.hstack((cellMassGrowthRateMultigen, cellmassGrowthRate))
			proteinGrowthRateMultigen = np.hstack((proteinGrowthRateMultigen, proteinGrowthRate))
			rnaGrowthRateMultigen = np.hstack((rnaGrowthRateMultigen, rnaGrowthRate))

		# seriesScrubber(cellMassGrowthRateMultigen, 1.5)
		# seriesScrubber(proteinGrowthRateMultigen, 1.5)
		# seriesScrubber(rnaGrowthRateMultigen, 1.5)


		cellMassGrowthRateMultigen = np.convolve(cellMassGrowthRateMultigen, np.ones(width) / width, mode = "same")
		proteinGrowthRateMultigen = np.convolve(proteinGrowthRateMultigen, np.ones(width) / width, mode = "same")
		rnaGrowthRateMultigen = np.convolve(rnaGrowthRateMultigen, np.ones(width) / width, mode = "same")
		# seriesScrubber(proteinGrowthRateMultigen, 2)

			# rnaGrowthRate = np.convolve(rnaGrowthRate, np.ones(width) / width, mode = "same")
			# dnaGrowthRate = np.convolve(dnaGrowthRate, np.ones(width) / width, mode = "same")

		# 	seriesScrubber(rnaGrowthRate,1.25)
			# seriesScrubber(dnaGrowthRate,3.25)

		linewidth = 2
		axis.plot(timeMultigen[:-width] / 60., cellMassGrowthRateMultigen[:-width] * 60., color = colors[0], alpha=0.9, label="Cell mass", linewidth=linewidth)
		axis.plot(timeMultigen[:-width] / 60., proteinGrowthRateMultigen[:-width] * 60., color = colors[1], alpha=0.9, label="Protein fraction", linewidth=linewidth)
		axis.plot(timeMultigen[:-width] / 60., rnaGrowthRateMultigen[:-width] * 60., color = colors[2], alpha=0.9, label="RNA fraction", linewidth=linewidth)
		axis.legend(loc=4,frameon=False)
		axis.set_ylim([0.014, 0.029])

		whitePadSparklineAxis(axis)
		axis.set_xlabel("Time (min)")
		axis.set_ylabel("Averaged instantaneous growth rate " + r"$(\frac{\frac{dX}{dt}}{X} [=] \frac{1}{min})$")
		plt.subplots_adjust(bottom = 0.2, left = 0.2)
		exportFigure(plt, plotOutDir, plotOutFileName,metadata)


		for axes in [axis]:
			axes.tick_params(
				axis='x',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				bottom=False,      # ticks along the bottom edge are off
				top=False,         # ticks along the top edge are off
				labelbottom=False) # labels along the bottom edge are off
			axes.tick_params(
				axis='y',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				left=False,        # ticks along the bottom edge are off
				right=False,       # ticks along the top edge are off
				labelleft=False)   # labels along the bottom edge are off

			axes.set_xlabel("")
			axes.set_ylabel("")

		plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped" ,metadata, transparent = True)
		plt.close("all")

		################ ALTERNATE SIZE ######################
		mult = 3

		fig, axis = plt.subplots(1,1)
		fig.set_figwidth(mm2inch(29.25)*mult)
		fig.set_figheight(mm2inch(24.401)*mult)

		width = 200

		timeMultigen = np.zeros(0)
		cellMassGrowthRateMultigen = np.zeros(0)
		proteinGrowthRateMultigen = np.zeros(0)
		rnaGrowthRateMultigen = np.zeros(0)


		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			time = main_reader.readColumn("time")
			timeStep = main_reader.readColumn("timeStepSec")

			mass = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = mass.readColumn("cellMass")
			proteinMass = mass.readColumn("proteinMass")
			rnaMass = mass.readColumn("rnaMass")
			dnaMass = mass.readColumn("dnaMass")

			cellmassGrowthRate = np.diff(cellMass) / cellMass[1:] / timeStep[:-1]
			proteinGrowthRate = np.diff(proteinMass) / proteinMass[1:] / timeStep[:-1]
			rnaGrowthRate = np.diff(rnaMass) / rnaMass[1:] / timeStep[:-1]
			dnaGrowthRate = np.diff(dnaMass) / dnaMass[1:] / timeStep[:-1]

			timeMultigen = np.hstack((timeMultigen, time[:-1]))
			cellMassGrowthRateMultigen = np.hstack((cellMassGrowthRateMultigen, cellmassGrowthRate))
			proteinGrowthRateMultigen = np.hstack((proteinGrowthRateMultigen, proteinGrowthRate))
			rnaGrowthRateMultigen = np.hstack((rnaGrowthRateMultigen, rnaGrowthRate))

		# seriesScrubber(cellMassGrowthRateMultigen, 1.5)
		# seriesScrubber(proteinGrowthRateMultigen, 1.5)
		# seriesScrubber(rnaGrowthRateMultigen, 1.5)


		cellMassGrowthRateMultigen = np.convolve(cellMassGrowthRateMultigen, np.ones(width) / width, mode = "same")
		proteinGrowthRateMultigen = np.convolve(proteinGrowthRateMultigen, np.ones(width) / width, mode = "same")
		rnaGrowthRateMultigen = np.convolve(rnaGrowthRateMultigen, np.ones(width) / width, mode = "same")
		# seriesScrubber(proteinGrowthRateMultigen, 2)

			# rnaGrowthRate = np.convolve(rnaGrowthRate, np.ones(width) / width, mode = "same")
			# dnaGrowthRate = np.convolve(dnaGrowthRate, np.ones(width) / width, mode = "same")

		# 	seriesScrubber(rnaGrowthRate,1.25)
			# seriesScrubber(dnaGrowthRate,3.25)

		linewidth = 2
		axis.plot(timeMultigen[:-width] / 60., cellMassGrowthRateMultigen[:-width] * 60., color = colors[0], alpha=0.9, label="Cell mass", linewidth=linewidth)
		axis.plot(timeMultigen[:-width] / 60., proteinGrowthRateMultigen[:-width] * 60., color = colors[1], alpha=0.9, label="Protein fraction", linewidth=linewidth)
		axis.plot(timeMultigen[:-width] / 60., rnaGrowthRateMultigen[:-width] * 60., color = colors[2], alpha=0.9, label="RNA fraction", linewidth=linewidth)
		# axis.legend(loc=4,frameon=False)
		axis.set_ylim([0.014, 0.029])

		whitePadSparklineAxis(axis)
		axis.set_xlabel("Time (min)")
		axis.set_ylabel("Averaged instantaneous growth rate " + r"$(\frac{\frac{dX}{dt}}{X} [=] \frac{1}{min})$")
		plt.subplots_adjust(bottom = 0.2, left = 0.2)
		exportFigure(plt, plotOutDir, plotOutFileName + "_altsize",metadata)

		for axes in [axis]:
			axes.tick_params(
				axis='x',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				bottom=False,      # ticks along the bottom edge are off
				top=False,         # ticks along the top edge are off
				labelbottom=False) # labels along the bottom edge are off
			axes.tick_params(
				axis='y',          # changes apply to the x-axis
				which='both',      # both major and minor ticks are affected
				left=False,        # ticks along the bottom edge are off
				right=False,       # ticks along the top edge are off
				labelleft=False)   # labels along the bottom edge are off

			axes.set_xlabel("")
			axes.set_ylabel("")

		plt.subplots_adjust(top = 1, bottom = trim_1, left = trim_1, right = 1)

		exportFigure(plt, plotOutDir, plotOutFileName + "_altsize_stripped" ,metadata, transparent = True)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
