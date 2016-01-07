#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir)

	# TODO: Declutter Y-axis

	# Get all cells
	firstCellLineage = []
	for gen_idx in range(ap.n_generations):
		firstCellLineage.append(ap.getGeneration(gen_idx)[0])

	fig, axesList = plt.subplots(2)

	## Create composition plot
	massNames = [
				#"dryMass",
				"proteinMass",
				"tRnaMass",
				"rRnaMass",
				'mRnaMass',
				"dnaMass",
				"smallMoleculeMass",
				]

	cleanNames = [
				#"Dry\nmass",
				"Protein",
				"tRNA",
				"rRNA",
				"mRNA",
				"DNA",
				"Small molecules"
				]

	colors = [
				"red",
				"blue",
				"green",
				"yellow",
				"grey",
				"cyan"
			]

	simDir = firstCellLineage[0] # Only considering the first generation for now
	time, massData = getMassData(simDir, massNames)
	massDataNorm = massData / massData.sum(axis = 0)

	initialComp = massDataNorm[:,0]
	finalComp = massDataNorm[:,-1]
	initFinalIdx = np.arange(2)
	barWidth = 0.75
	oldLayerData = np.zeros(2)

	for idx, massType in enumerate(massNames):
		layerData = np.array([initialComp[idx], finalComp[idx]])
		axesList[0].bar(initFinalIdx, layerData, barWidth, bottom=oldLayerData, label=cleanNames[idx], color=colors[idx])
		oldLayerData += layerData
	axesList[0].legend(prop={'size':6})
	axesList[0].set_title('Mass composition')
	axesList[0].set_xticks(initFinalIdx + barWidth/2.)
	axesList[0].set_xticklabels(("$t_i=${}".format(time[0]), "$t_f=${}".format(time[-1])))

	## Create log growth plot
	massNames = [
				"dryMass",
				"proteinMass",
				"rnaMass",
				"dnaMass",
				]

	cleanNames = [
				"Dry",
				"Protein",
				"RNA",
				"DNA",
				]

	colors = [
				"black",
				"red",
				"blue",
				"green",
			]
	simDir = firstCellLineage[0] # Only considering the first generation for now
	time, massData = getMassData(simDir, massNames)
	for idx, massType in enumerate(massNames):
		axesList[1].plot(time / 60., np.log(massData[idx,:]), label=cleanNames[idx], color=colors[idx])
	axesList[1].legend(prop={'size':6})
	axesList[1].set_xlabel('Time (min)')
	axesList[1].set_ylabel('log(mass in fg)')

	# import ipdb; ipdb.set_trace()
	# for simDir in firstCellLineage:
	# 	simOutDir = os.path.join(simDir, "simOut")
	# 	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	# 	mass = TableReader(os.path.join(simOutDir, "Mass"))

	# 	massData = np.zeros((len(massNames),time.size))

	# 	for idx, massType in enumerate(massNames):
	# 		massData[idx,:] = mass.readColumn(massNames[idx])

	# 	massData = massData / massData.sum(axis = 0)

	# 	for idx, massType in enumerate(massNames):
	# 		axesList[idx].plot(time / 60, massData[idx,:])
	# 		axesList[idx].set_ylabel(cleanNames[idx])

	# for axes in axesList:
	# 	axes.set_yticks(list(axes.get_ylim()))

	# axesList[-1].set_xlabel('Time (min)')

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
	plt.close("all")


def getMassData(simDir, massNames):
	simOutDir = os.path.join(simDir, "simOut")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	mass = TableReader(os.path.join(simOutDir, "Mass"))

	massData = np.zeros((len(massNames),time.size))

	for idx, massType in enumerate(massNames):
		massData[idx,:] = mass.readColumn(massNames[idx])

	return time, massData


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
