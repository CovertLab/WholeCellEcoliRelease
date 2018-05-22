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

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"
	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	sim_data = cPickle.load(open(simDataFile, "rb"))

	T_ADD_AA = None
	T_CUT_AA = None
	nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
	if "aa" in nutrientsTimeSeriesLabel:
		if "add" in nutrientsTimeSeriesLabel and "cut" in nutrientsTimeSeriesLabel:
			T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0]
			T_CUT_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][2][0]

	# Get all cells
	allDir = ap.get_cells()
	nCells = allDir.shape[0]
	nGens = ap.n_generation

	massNames = ["dryMass", "proteinMass", "rnaMass", "dnaMass",]
	cleanNames = ["Dry mass", "Protein mass", "RNA Mass", "DNA mass",]

	fig = plt.figure(figsize = (14, 10))
	ax1 = plt.subplot2grid((4,2), (0,0), rowspan = 4)
	ax2 = plt.subplot2grid((4,2), (0,1))
	ax3 = plt.subplot2grid((4,2), (1,1))
	ax4 = plt.subplot2grid((4,2), (2,1))
	ax5 = plt.subplot2grid((4,2), (3,1))
	axesList = [ax2, ax3, ax4, ax5]
	colors = ["blue", "green", "red", "cyan"]

	for gen, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx])
			ax1.semilogy(time / 60. / 60., (2. ** gen) * massToPlot, linewidth = 2, color = colors[idx])
			dMassToPlot = (np.log10((2. ** gen) * massToPlot[1:]) - np.log10((2. ** gen) * massToPlot[:-1]))
			axesList[idx].semilogy(time[1 + NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE:-NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE] / 60. / 60., dMassToPlot[NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE:-NUM_SKIP_TIMESTEPS_AT_GEN_CHANGE], color = colors[idx])

			# Todo: plot exponential fit to intial condition

	ax1.legend(cleanNames, loc = "best", fontsize = 10)
	ax1.set_ylabel("Mass (fg)")
	ax1.set_xlabel("Time (hour)")
	if T_ADD_AA != None and T_CUT_AA !=None:
		ax1.plot((T_ADD_AA / 60. / 60., T_ADD_AA / 60. / 60.), ax1.get_ylim(), "k--")
		ax1.plot((T_CUT_AA / 60. / 60., T_CUT_AA / 60. / 60.), ax1.get_ylim(), "k--")

	ax2.set_title("Instantaneous change in mass (fg / second)")
	ax2.set_xlim((-0.1, ax2.get_xlim()[1]))
	ax5.set_xlabel("Time (hour)")
	for idx, ax in enumerate(axesList):
		ax.set_ylabel(cleanNames[idx])
		ax.set_xlim(0, time[-1] / 60. / 60.)
		if T_ADD_AA != None and T_CUT_AA !=None:
			ax.plot((T_ADD_AA / 60. / 60., T_ADD_AA / 60. / 60.), ax.get_ylim(), "k--")
			ax.plot((T_CUT_AA / 60. / 60., T_CUT_AA / 60. / 60.), ax.get_ylim(), "k--")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.2)

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
