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
from matplotlib import pyplot as plt
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import curve_fit
import wholecell.utils.constants

NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE = 1

from wholecell.utils.sparkline import whitePadSparklineAxis

trim = 0.03
trim_1 = 0.06
def mm2inch(value):
	return value * 0.0393701

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
	nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label
	if "aa" in nutrients_time_series_label:
		if "add" in nutrients_time_series_label and "cut" in nutrients_time_series_label:
			T_ADD_AA = sim_data.external_state.environment.nutrients_time_series[nutrients_time_series_label][1][0]
			T_CUT_AA = sim_data.external_state.environment.nutrients_time_series[nutrients_time_series_label][2][0]

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
	axis.set_ylim([0.014, 0.029])

	whitePadSparklineAxis(axis)
	axis.set_xlabel("Time (min)")
	axis.set_ylabel("Averaged instantaneous growth rate " + r"$(\frac{\frac{dX}{dt}}{X} [=] \frac{1}{min})$")
	plt.subplots_adjust(bottom = 0.2, left = 0.2)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)


	for axes in [axis]:
		axes.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
		axes.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			left='off',      # ticks along the bottom edge are off
			right='off',         # ticks along the top edge are off
			labelleft='off') # labels along the bottom edge are off

		axes.set_xlabel("")
		axes.set_ylabel("")

	plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

	from wholecell.analysis.analysis_tools import exportFigure
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
	# axis.legend(loc=4,frameon=False)
	axis.set_ylim([0.014, 0.029])

	whitePadSparklineAxis(axis)
	axis.set_xlabel("Time (min)")
	axis.set_ylabel("Averaged instantaneous growth rate " + r"$(\frac{\frac{dX}{dt}}{X} [=] \frac{1}{min})$")
	plt.subplots_adjust(bottom = 0.2, left = 0.2)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_altsize",metadata)


	for axes in [axis]:
		axes.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
		axes.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			left='off',      # ticks along the bottom edge are off
			right='off',         # ticks along the top edge are off
			labelleft='off') # labels along the bottom edge are off

		axes.set_xlabel("")
		axes.set_ylabel("")

	plt.subplots_adjust(top = 1, bottom = trim_1, left = trim_1, right = 1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_altsize_stripped" ,metadata, transparent = True)
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
