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
import cPickle

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# Get first cell from each generation
	firstCellLineage = []
	for gen_idx in range(ap.n_generation):
		firstCellLineage.append(ap.get_cells(generation = [gen_idx])[0])

	sim_data = cPickle.load(open(simDataFile, "rb"))
	max_elongationRate = 21. # TODO: Fix this
	elongationRate = float(sim_data.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s))

	fig = plt.figure()
	fig.set_size_inches(10,12)

	## Create composition plot
	massNames = [
				"proteinMass",
				"rnaMass",
				]

	cleanNames = [
				"Protein",
				"RNA",
				]

	colors = [
				"red",
				"blue",
			]

	# Composition plotting
	# Create plotting indexes with some offsets to make it look nice
	initFinalIdx = range(4 * len(firstCellLineage))
	for idx, value in enumerate(initFinalIdx):
		if idx % 2:
			initFinalIdx[idx] = initFinalIdx[idx] - 0.1
		else:
			initFinalIdx[idx] = initFinalIdx[idx] + 0.1
	initFinalIdx = np.array(initFinalIdx)
	distanceInitFinalIdx = np.arange(start = 0.5, stop = 4 * len(firstCellLineage) + 0.5, step = 2)
	tickLabels = tuple(["$t=0$\nActual","Expected","$t=t_d$\nActual","Expected"] * len(firstCellLineage))
	barWidth = 0.75
	gs = gridspec.GridSpec(6, 4)
	ax1 = plt.subplot(gs[:1, :-2])
	ax8 = plt.subplot(gs[1, :-2])
	ax2 = plt.subplot(gs[:2, 2:])


	for gen, simDir in enumerate(firstCellLineage):
		simOutDir = os.path.join(simDir, "simOut")

		# Get simulation composition
		time, massData = getMassData(simDir, massNames)
		timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
		_, dryMass = getMassData(simDir, ["dryMass"])

		# massDataNorm = massData / massData.sum(axis = 0)
		massDataNorm = massData / dryMass

		initialComp = massDataNorm[:,1:6*60+1].mean(axis=1)
		finalComp = massDataNorm[:,-5*60:].mean(axis=1)
		initialDryMass = dryMass[0]
		finalDryMass = dryMass[-1]

		# Get expected composition based on growth rates
		_, growthRate = getMassData(simDir, ["instantaniousGrowthRate"])

		initialGrowthRate = np.mean(growthRate[1:6*60+1]) * 60
		initialDoublingTime = np.log(2) / initialGrowthRate * units.min # Five mintues skipping the first minute
		finalGrowthRate = (np.mean(growthRate[-5*60:]) * 60)
		finalDoublingTime = np.log(2) / finalGrowthRate * units.min # Last five minutes

		expectedInitialComp, expectedInitialMass = getExpectedComposition(initialDoublingTime)
		#expectedInitialComp = expectedInitialComp / expectedInitialComp.sum()
		expectedInitialComp = expectedInitialComp / expectedInitialMass.asNumber(units.fg)
		expectedFinalComp, finalCompInitMass = getExpectedComposition(finalDoublingTime)
		#expectedFinalComp = expectedFinalComp / expectedFinalComp.sum()
		expectedFinalMass = finalCompInitMass * 2
		expectedFinalComp = expectedFinalComp / expectedFinalMass.asNumber(units.fg) * 2
		initialDistance = np.linalg.norm(expectedInitialComp - initialComp, 2) / np.linalg.norm(expectedInitialComp, 2)
		finalDistance = np.linalg.norm(expectedFinalComp - finalComp, 2) / np.linalg.norm(expectedFinalComp, 2)

		# Calculate ribosomal rna doubling times
		rrn16S_produced = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn16S_produced")
		rrn23S_produced = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn23S_produced")
		rrn5S_produced = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn5S_produced")

		ids_16s = []
		ids_16s.extend(sim_data.moleculeGroups.s30_16sRRNA)
		ids_16s.extend(sim_data.moleculeGroups.s30_fullComplex)

		ids_23s = []
		ids_23s.extend(sim_data.moleculeGroups.s50_23sRRNA)
		ids_23s.extend(sim_data.moleculeGroups.s50_fullComplex)

		ids_5s = []
		ids_5s.extend(sim_data.moleculeGroups.s50_5sRRNA)
		ids_5s.extend(sim_data.moleculeGroups.s50_fullComplex)
		
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")

		idx_16s = np.array([moleculeIds.index(comp) for comp in ids_16s], np.int)
		idx_23s = np.array([moleculeIds.index(comp) for comp in ids_23s], np.int)
		idx_5s = np.array([moleculeIds.index(comp) for comp in ids_5s], np.int)

		rrn16s_count_bulk = bulkMolecules.readColumn("counts")[:, idx_16s].sum(axis=1)
		rrn23s_count_bulk = bulkMolecules.readColumn("counts")[:, idx_23s].sum(axis=1)
		rrn5s_count_bulk = bulkMolecules.readColumn("counts")[:, idx_5s].sum(axis=1)

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		ribosome_count_unique = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		rrn16s_count = rrn16s_count_bulk + ribosome_count_unique
		rrn23s_count = rrn23s_count_bulk + ribosome_count_unique
		rrn5s_count = rrn5s_count_bulk + ribosome_count_unique

		rrn16S_doubling_time = units.s * np.log(2) / (rrn16S_produced / timeStep / rrn16s_count)
		rrn23S_doubling_time = units.s * np.log(2) / (rrn23S_produced / timeStep / rrn23s_count)
		rrn5S_doubling_time = units.s * np.log(2) / (rrn5S_produced / timeStep / rrn5s_count)

		rrn16S_doubling_time[rrn16S_doubling_time.asNumber() == np.inf] = np.nan * units.s
		rrn23S_doubling_time[rrn23S_doubling_time.asNumber() == np.inf] = np.nan * units.s
		rrn5S_doubling_time[rrn5S_doubling_time.asNumber() == np.inf] = np.nan * units.s

		# Plotting

		oldLayerData = np.zeros(4)
		indexesForGen = initFinalIdx[gen*4:(gen+1)*4]
		ax1.axvline(x = indexesForGen[-1] + barWidth + 0.225, linewidth=2, color='k', linestyle='--')

		for idx, massType in enumerate(massNames):
			layerData = np.array([initialComp[idx], expectedInitialComp[idx], finalComp[idx], expectedFinalComp[idx]])
			if gen == 0:
				ax1.bar(left=indexesForGen, height=layerData, width=barWidth, bottom=oldLayerData, label=cleanNames[idx], color=colors[idx])
			else:
				ax1.bar(left=indexesForGen, height=layerData, width=barWidth, bottom=oldLayerData, color=colors[idx])
			oldLayerData += layerData
		ax1.set_ylim([0, 1])
		ax1.legend(prop={'size':4})
		ax1.set_title('Mass composition')

		# Plot Euclidian distance
		distanceIdxForGen = distanceInitFinalIdx[gen*2:(gen+1)*2]
		ax8.bar(distanceIdxForGen, [initialDistance, finalDistance])
		ax8.set_ylabel("Relative error")

		# Plot dry mass
		ax2.bar(indexesForGen, [initialDryMass, expectedInitialMass.asNumber(units.fg), finalDryMass, expectedFinalMass.asNumber(units.fg)], barWidth)
		ax2.set_title('Dry mass (fg)')
		ax2.axvline(x = indexesForGen[-1] + barWidth + 0.225, linewidth=2, color='k', linestyle='--')

		# Plot growth rate
		ax3 = plt.subplot(gs[2, :-1])
		doublingTime = np.log(2) / growthRate / 60. # min
		avgDoublingTime = doublingTime[1:].mean()
		stdDoublingTime = doublingTime[1:].std()
		ax3.plot(time / 60., doublingTime)
		ax3.set_ylim([0, avgDoublingTime + 2*stdDoublingTime])
		ax3.set_ylabel("Doubling time (min)")
		ax3.axvline(x = time.max() / 60., linewidth=2, color='k', linestyle='--')

		# Plot ribosome data
		simOutDir = os.path.join(simDir, "simOut")
		ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
		effectiveElongationRate = ribosomeDataFile.readColumn("effectiveElongationRate")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		ribosomeDataFile.close()

		ax4 = plt.subplot(gs[3, :-1])
		ax4.plot(time / 60., effectiveElongationRate, label="Effective elongation rate", linewidth=2)
		ax4.plot(time / 60., max_elongationRate * np.ones(time.size), 'r--')
		ax4.set_ylabel("Effective elongation\nrate (aa/s/ribosome)")
		ax4.axvline(x = time.max() / 60., linewidth=2, color='k', linestyle='--')

		# Plot rrn counts
		# bulkMoleculesFile = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		# bulkMoleculeIds = bulkMoleculesFile.readAttribute("objectNames")
		# rrn_idx = bulkMoleculesFile.readAttribute("objectNames").index('rrn_operon')
		# rrn_counts = bulkMoleculesFile.readColumn("counts")[:, rrn_idx]
		# bulkMoleculesFile.close()

		# ax5 = plt.subplot(gs[4, :-1])
		# ax5.plot(time / 60., rrn_counts, label="Rrn operon counts", linewidth=2)
		# ax5.set_ylim([rrn_counts.min() - 1, rrn_counts.max() + 1])
		# ax5.set_ylabel("Rrn operons")
		# ax5.axvline(x = time.max() / 60., linewidth=2, color='k', linestyle='--')

		ax6 = plt.subplot(gs[5, :-1])
		ax6.plot(time / 60., rrn16S_doubling_time.asNumber(units.min), linewidth=2)
		ax6.plot(time / 60., rrn23S_doubling_time.asNumber(units.min), linewidth=2)
		ax6.plot(time / 60., rrn5S_doubling_time.asNumber(units.min), linewidth=2)
		ax6.axvline(x = time.max() / 60., linewidth=2, color='k', linestyle='--')
		ax6.set_ylabel("rrn doubling time")
		ax6.set_xlabel("Time (min)")

	# Set axes labels for x-axis in plots 1 and 2
	ax8.set_xticks(initFinalIdx + barWidth/2.)
	ax8.set_xticklabels(tickLabels)
	plt.setp(ax8.xaxis.get_majorticklabels(), rotation = 70, fontsize=8)
	ax2.set_xticks(initFinalIdx + barWidth/2.)
	ax2.set_xticklabels(tickLabels)
	plt.setp(ax2.xaxis.get_majorticklabels(), rotation = 70, fontsize=8)

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



	ax7 = plt.subplot(gs[2:, -1])
	simDir = firstCellLineage[0] # Only considering the first generation for now
	time, massData = getMassData(simDir, massNames)
	for idx, massType in enumerate(massNames):
		ax7.plot(time / 60., np.log(massData[idx,:]), label=cleanNames[idx], color=colors[idx])
	ax7.legend(prop={'size':6})
	ax7.set_xlabel('Time (min)')
	ax7.set_ylabel('log(mass in fg)')

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
	if doubling_time < 24. * units.min or doubling_time > 100. * units.min:
		return np.zeros(2), 0 * units.fg

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

	masses = np.array([protein, rRnaMass + tRNA + mRnaMass]) / sim_data.mass.avgCellToInitialCellConvFactor
	
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
