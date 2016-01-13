#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir)

	# Get all cells
	firstCellLineage = []
	for gen_idx in range(ap.n_generations):
		firstCellLineage.append(ap.getGeneration(gen_idx)[0])

	# fig, axesList = plt.subplots(3)
	# fig.set_size_inches(10,20)

	fig = plt.figure()
	fig.set_size_inches(10,10)

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

	# Get simulation composition
	simDir = firstCellLineage[0] # Only considering the first generation for now
	time, massData = getMassData(simDir, massNames)
	_, dryMass = getMassData(simDir, ["dryMass"])

	massDataNorm = massData / massData.sum(axis = 0)

	initialComp = massDataNorm[:,0]
	finalComp = massDataNorm[:,-1]
	initialDryMass = dryMass[0]
	finalDryMass = dryMass[-1]

	# Get expected composition based on growth rates
	_, growthRate = getMassData(simDir, ["instantaniousGrowthRate"])

	if time.size < 6 * 60:
		print "Simulation isn't long enough, exiting."
		return

	initialGrowthRate = np.mean(growthRate[1*60:6*60]) * 60
	initialDoublingTime = np.log(2) / initialGrowthRate * units.min # Five mintues skipping the first minute
	finalGrowthRate = (np.mean(growthRate[-5*60:]) * 60)
	finalDoublingTime = np.log(2) / finalGrowthRate * units.min # Last five minutes

	expectedInitialComp, expectedInitialMass = getExpectedComposition(initialDoublingTime)
	expectedInitialComp = expectedInitialComp / expectedInitialComp.sum()
	expectedFinalComp, finalCompInitMass = getExpectedComposition(finalDoublingTime)
	expectedFinalComp = expectedFinalComp / expectedFinalComp.sum()
	expectedFinalMass = finalCompInitMass * 2

	initialDistance = np.linalg.norm(expectedInitialComp - initialComp, 2)
	finalDistance = np.linalg.norm(expectedFinalComp - finalComp, 2)
	distanceInitialToFinal = np.linalg.norm(initialComp - finalComp, 2)

	initFinalIdx = np.arange(4)
	barWidth = 0.75
	oldLayerData = np.zeros(4)

	gs = gridspec.GridSpec(4, 4)

	ax1 = plt.subplot(gs[:-2, :-2])

	for idx, massType in enumerate(massNames):
		layerData = np.array([initialComp[idx], expectedInitialComp[idx], finalComp[idx], expectedFinalComp[idx]])
		ax1.bar(initFinalIdx, layerData, barWidth, bottom=oldLayerData, label=cleanNames[idx], color=colors[idx])
		oldLayerData += layerData
	ax1.legend(prop={'size':6})
	ax1.set_title('Mass composition')
	ax1.set_xticks(initFinalIdx + barWidth/2.)
	ax1.set_xticklabels(("{}".format(time[0]), "E({})".format(time[0]), "{}".format(time[-1]), "E({})".format(time[-1])))
	ax1.set_xlabel("Time (s)")

	ax2 = plt.subplot(gs[:-2, 2])
	ax2.bar(initFinalIdx, [initialDryMass, expectedInitialMass.asNumber(units.fg), finalDryMass, expectedFinalMass.asNumber(units.fg)], barWidth)
	ax2.set_title('Dry mass (fg)')
	ax2.set_xticks(initFinalIdx + barWidth/2.)
	ax2.set_xticklabels(("{}".format(time[0]), "E({})".format(time[0]), "{}".format(time[-1]), "E({})".format(time[-1])))
	ax2.set_xlabel("Time (s)")

	ax3 = plt.subplot(gs[:-2, -1])
	ax3.bar(np.arange(2), [initialGrowthRate * 60, finalGrowthRate * 60], barWidth)
	ax3.set_title('Growth rate (1/hr)')
	ax3.set_xticks(np.arange(2) + barWidth/2.)
	ax3.set_xticklabels(("{}".format(time[0]), "{}".format(time[-1])))
	ax3.set_xlabel("Time (s)")

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

	ax4 = plt.subplot(gs[2:, -1])
	simDir = firstCellLineage[0] # Only considering the first generation for now
	time, massData = getMassData(simDir, massNames)
	for idx, massType in enumerate(massNames):
		ax4.plot(time / 60., np.log(massData[idx,:]), label=cleanNames[idx], color=colors[idx])
	ax4.legend(prop={'size':6})
	ax4.set_xlabel('Time (min)')
	ax4.set_ylabel('log(mass in fg)')

	ax5 = plt.subplot(gs[2, :-1])
	ax5.plot(time / 60., growthRate * 3600)
	ax5.set_ylim([0 * 3600, 0.0003 * 3600])
	ax5.set_title("Growth rate (1/hr)")
	ax5.set_xlabel("Time (min)")

	fig.subplots_adjust(hspace=.5, wspace = 0.3)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
	plt.close("all")


def getMassData(simDir, massNames):
	simOutDir = os.path.join(simDir, "simOut")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	mass = TableReader(os.path.join(simOutDir, "Mass"))

	massFractionData = np.zeros((len(massNames),time.size))

	for idx, massType in enumerate(massNames):
		massFractionData[idx,:] = mass.readColumn(massNames[idx])

	if len(massNames) == 1:
		massFractionData = massFractionData.reshape(-1)

	return time, massFractionData

def getExpectedComposition(doubling_time):
	#return np.ones(6), 100. * units.fg

	from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
	raw_data = KnowledgeBaseEcoli()

	from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
	sim_data = fitSimData_1(raw_data, doubling_time = doubling_time)

	subMasses = sim_data.mass.avgCellSubMass
	protein = subMasses['proteinMass'].asNumber(units.fg)
	tRNA = subMasses['tRnaMass'].asNumber(units.fg)
	rRnaMass = (subMasses['rRna23SMass'] + subMasses['rRna5SMass'] + subMasses['rRna16SMass']).asNumber(units.fg)
	mRnaMass = subMasses['mRnaMass'].asNumber(units.fg)
	dnaMass = subMasses['dnaMass'].asNumber(units.fg)
	smallMolecules = sim_data.mass.fitAvgSolublePoolMass.asNumber(units.fg)

	masses = np.array([protein, tRNA, rRnaMass, mRnaMass, dnaMass, smallMolecules]) / sim_data.mass.avgCellToInitialCellConvFactor

	initialMass = sim_data.mass.avgCellDryMassInit

	return masses, initialMass

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
