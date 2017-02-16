#!/usr/bin/env python
"""
@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/19/2017
"""

import argparse
import os
import numpy as np
import cPickle
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import curve_fit
import wholecell.utils.constants

NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE = 1

from wholecell.utils.sparkline import whitePadSparklineAxis


def seriesScrubber(series, factor):
	series[abs(series - np.median(series)) > factor * np.nanstd(series)] = np.nan


def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"
	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, cohort_plot = True)
	sim_data = cPickle.load(open(simDataFile, "rb"))

	T_ADD_AA = None
	T_CUT_AA = None
	nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
	if "aa" in nutrientsTimeSeriesLabel:
		if "add" in nutrientsTimeSeriesLabel and "cut" in nutrientsTimeSeriesLabel:
			T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0]
			T_CUT_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][2][0]

	# Get all cells
	allDir = ap.get_cells(seed=[0])
	nCells = allDir.shape[0]
	nGens = ap.n_generation

	massNames = ["dryMass", "proteinMass", "rnaMass", "dnaMass",]
	cleanNames = ["Dry mass", "Protein mass", "RNA Mass", "DNA mass",]

	# fig = plt.figure(figsize = (14, 10))
	# ax1 = plt.subplot2grid((4,2), (0,0), rowspan = 4)
	# ax2 = plt.subplot2grid((4,2), (0,1))
	# ax3 = plt.subplot2grid((4,2), (1,1))
	# ax4 = plt.subplot2grid((4,2), (2,1))
	# ax5 = plt.subplot2grid((4,2), (3,1))
	# axesList = [ax2, ax3, ax4, ax5]
	# colors = ["blue", "green", "red", "cyan"]
	colors = ["#43aa98", "#0071bb", "#bf673c"]

	fig, axis = plt.subplots(1,1)
	fig.set_figwidth(10)
	fig.set_figheight(10)

	width = 200

	timeMultigen = np.zeros(0)
	cellMassGrowthRateMultigen = np.zeros(0)
	proteinGrowthRateMultigen = np.zeros(0)
	rnaGrowthRateMultigen = np.zeros(0)


	for gen, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
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
	# axis.set_ylim([0.00020, 0.00035])

	whitePadSparklineAxis(axis)
	axis.set_xlabel("Time (min)")
	axis.set_ylabel("Averaged instantanious growth rate " + r"$(\frac{\frac{dX}{dt}}{X} [=] \frac{1}{min})$")
	plt.subplots_adjust(bottom = 0.2, left = 0.2)
		# axis.plot(time[:-1] / 60., rnaGrowthRate, color = colors[1], alpha=0.5)
		# axis.plot(time[:-1] / 60., dnaGrowthRate, color = colors[2], alpha=0.5)

		# for idx, massType in enumerate(massNames):
		# 	massToPlot = mass.readColumn(massNames[idx])
		# 	dMassToPlot = ((2. ** gen) * massToPlot[1:] - (2. ** gen) * massToPlot[:-1])
		# 	axis.plot(time[1 + NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE:-NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE] / 60. / 60., dMassToPlot[NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE:-NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE], color = colors[idx])

		# 	# Todo: plot exponential fit to intial condition

	# ax1.legend(cleanNames, loc = "best", fontsize = 10)
	# ax1.set_ylabel("Mass (fg)")
	# ax1.set_xlabel("Time (hour)")
	# if T_ADD_AA != None and T_CUT_AA !=None:
	# 	ax1.plot((T_ADD_AA / 60. / 60., T_ADD_AA / 60. / 60.), ax1.get_ylim(), "k--")
	# 	ax1.plot((T_CUT_AA / 60. / 60., T_CUT_AA / 60. / 60.), ax1.get_ylim(), "k--")

	# ax2.set_title("Instantaneous change in mass (fg / second)")
	# ax2.set_xlim((-0.1, ax2.get_xlim()[1]))
	# ax5.set_xlabel("Time (hour)")
	# for idx, ax in enumerate(axesList):
	# 	ax.set_ylabel(cleanNames[idx])
	# 	ax.set_xlim(0, time[-1] / 60. / 60.)
	# 	if T_ADD_AA != None and T_CUT_AA !=None:
	# 		ax.plot((T_ADD_AA / 60. / 60., T_ADD_AA / 60. / 60.), ax.get_ylim(), "k--")
	# 		ax.plot((T_CUT_AA / 60. / 60., T_CUT_AA / 60. / 60.), ax.get_ylim(), "k--")

	# plt.subplots_adjust(hspace = 0.2, wspace = 0.2)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])