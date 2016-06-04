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

FONT_SIZE = 9
DOWN_SAMPLE = 100

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	print "Disabled. Calls fitter with doubling_time argument which is deprecated."
	return
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(variantDir, cohort_plot = True)
	allCells = ap.get_cells()

	fig = plt.figure()
	fig.set_size_inches(10,12)

	# Composition plotting
	gs = gridspec.GridSpec(4, 3)
	doublingTime_axis = plt.subplot(gs[0,0])
	rnaFrac_axis = plt.subplot(gs[1,0])
	proteinFrac_axis = plt.subplot(gs[2,0])
	dryMassInit_axis = plt.subplot(gs[3,0])
	dryMassFinal_axis = plt.subplot(gs[3,1])

	doublingTimeViolin_axis = plt.subplot(gs[0, 1:])
	rnaFractionViolin_axis = plt.subplot(gs[1, 1:])
	proteinFractionViolin_axis = plt.subplot(gs[2, 1:])


	sim_data = cPickle.load(open(simDataFile, "rb"))
	expectedProtein, expectedRna, expectedDryMassInit = getExpectedComposition(sim_data.doubling_time)

	maxTime = getMaxTime(allCells)

	gigsPerArray = maxTime * allCells.size * 8 / 1e9

	## Load values ##
	time = getSingleValue(allCells, tableName = "Main", colName = "time", maxTime = maxTime)
	growthRate = getSingleValue(allCells, tableName = "Mass", colName = "instantaniousGrowthRate", maxTime = maxTime)
	rnaMass = getSingleValue(allCells, tableName = "Mass", colName = "rnaMass", maxTime = maxTime)
	proteinMass = getSingleValue(allCells, tableName = "Mass", colName = "proteinMass", maxTime = maxTime)
	dryMass = getSingleValue(allCells, tableName = "Mass", colName = "dryMass", maxTime = maxTime)

	## Create growth rate histogram ##
	doublingTime = units.s * np.log(2) / growthRate

	histDoublingTime = removeNanReshape(doublingTime.asNumber(units.min))

	nbins = np.ceil(np.sqrt(histDoublingTime.size))
	doublingTime_axis.hist(histDoublingTime, nbins)
	doublingTime_axis.set_xlim([histDoublingTime.mean() - 3*histDoublingTime.std(), histDoublingTime.mean() + 3*histDoublingTime.std()])
	doublingTime_axis.axvline(x = sim_data.doubling_time.asNumber(units.min), linewidth=2, color='k', linestyle='--')
	doublingTime_axis.set_ylabel("Frequency of time point\nwith doubling time", fontsize=FONT_SIZE)
	doublingTime_axis.set_title("Instantanious doubling time (min) for\n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Create rna fraction histogram ##
	rnaFraction = rnaMass / dryMass
	expectedRnaFraction = expectedRna / expectedDryMassInit

	histRnaFraction = removeNanReshape(rnaFraction)

	nbins = np.ceil(np.sqrt(histRnaFraction.size))
	rnaFrac_axis.hist(histRnaFraction, nbins)
	features = np.array([histRnaFraction.mean() - 3*histRnaFraction.std(), expectedRnaFraction, histRnaFraction.mean() + 3*histRnaFraction.std()])
	x_min = features.min()
	x_max = features.max()
	rnaFrac_axis.set_xlim([x_min, x_max])
	rnaFrac_axis.set_xlim([histRnaFraction.mean() - 3*histRnaFraction.std(), histRnaFraction.mean() + 3*histRnaFraction.std()])
	rnaFrac_axis.axvline(x = expectedRnaFraction, linewidth=2, color='k', linestyle='--')
	rnaFrac_axis.set_ylabel("Frequency of time point\nwith fraction", fontsize=FONT_SIZE)
	rnaFrac_axis.set_title("Instantanious rna dry mass fraction for\n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Create protein fraction histogram ##
	proteinFraction = proteinMass / dryMass
	expectedProteinFraction = expectedProtein / expectedDryMassInit

	histProteinFraction = removeNanReshape(proteinFraction)

	nbins = np.ceil(np.sqrt(histProteinFraction.size))
	proteinFrac_axis.hist(histProteinFraction, nbins)
	features = np.array([histProteinFraction.mean() - 3*histProteinFraction.std(), expectedProteinFraction, histProteinFraction.mean() + 3*histProteinFraction.std()])
	x_min = features.min()
	x_max = features.max()
	proteinFrac_axis.set_xlim([x_min, x_max])
	proteinFrac_axis.axvline(x = expectedProteinFraction, linewidth=2, color='k', linestyle='--')
	proteinFrac_axis.set_ylabel("Frequency of time point\nwith fraction", fontsize=FONT_SIZE)
	proteinFrac_axis.set_title("Instantanious protein dry mass fraction for\n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Create dry mass initial and final plots ##
	dryMassInit = dryMass[:,0]
	dryMassFinal = np.array([dryMass[idx,:][~np.isnan(dryMass[idx,:])][-1] for idx in range(dryMass.shape[0])])
	expectedDryMassFinal = 2 * expectedDryMassInit
	
	nbins = np.ceil(np.sqrt(dryMassInit.size))
	dryMassInit_axis.hist(dryMassInit, nbins)
	features = np.array([dryMassInit.mean() - 3*dryMassInit.std(), expectedDryMassInit, dryMassInit.mean() + 3*dryMassInit.std()])
	x_min = features.min()
	x_max = features.max()
	dryMassInit_axis.set_xlim([x_min, x_max])
	dryMassInit_axis.axvline(x = expectedDryMassInit, linewidth=2, color='k', linestyle='--')
	dryMassInit_axis.set_ylabel("Frequency of cell\nwith dry mass", fontsize=FONT_SIZE)
	dryMassInit_axis.set_title("Initial dry mass (fg) for \n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	nbins = np.ceil(np.sqrt(dryMassFinal.size))
	dryMassFinal_axis.hist(dryMassFinal, nbins)
	features = np.array([dryMassFinal.mean() - 3*dryMassFinal.std(), expectedDryMassFinal, dryMassFinal.mean() + 3*dryMassFinal.std()])
	x_min = features.min()
	x_max = features.max()
	dryMassFinal_axis.set_xlim([x_min, x_max])
	dryMassFinal_axis.axvline(x = expectedDryMassFinal, linewidth=2, color='k', linestyle='--')
	dryMassFinal_axis.set_ylabel("Frequency of cell\nwith dry mass", fontsize=FONT_SIZE)
	dryMassFinal_axis.set_title("Final dry mass (fg) for \n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Create doubling time violin plots ##
	downSampleDoublingTime = downSample(doublingTime.asNumber(units.min), DOWN_SAMPLE)
	positions = np.arange(start = 0, stop = doublingTime.asNumber().shape[1], step = np.floor(doublingTime.asNumber().shape[1] / downSampleDoublingTime.shape[1]))[:downSampleDoublingTime.shape[1]]
	width = positions.max() / positions.size
	doublingTimeViolin_axis.violinplot(downSampleDoublingTime, widths=width, showmeans=True, positions = positions)
	doublingTimeViolin_axis.set_title("Instantanious doubling time (min) for\n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Create rna fraction violin plots ##
	downSampleRnaFraction = downSample(rnaFraction, DOWN_SAMPLE)
	positions = np.arange(start = 0, stop = rnaFraction.shape[1], step = np.floor(rnaFraction.shape[1] / downSampleRnaFraction.shape[1]))[:downSampleRnaFraction.shape[1]]
	width = positions.max() / positions.size
	rnaFractionViolin_axis.violinplot(downSampleRnaFraction, widths=width, showmeans=True, positions = positions)
	rnaFractionViolin_axis.set_title("Rna dry mass fraction for\n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Create protein fraction violin plots ##
	downSampleProteinFraction = downSample(proteinFraction, DOWN_SAMPLE)
	positions = np.arange(start = 0, stop = proteinFraction.shape[1], step = np.floor(proteinFraction.shape[1] / downSampleProteinFraction.shape[1]))[:downSampleProteinFraction.shape[1]]
	width = positions.max() / positions.size
	proteinFractionViolin_axis.violinplot(downSampleProteinFraction, widths=width, showmeans=True, positions = positions)
	proteinFractionViolin_axis.set_title("Protein dry mass fraction for\n{} generations of cells (n={})".format(ap.n_generation, len(allCells)), fontsize=FONT_SIZE)

	## Set final formatting

	fig.subplots_adjust(hspace=.3, wspace = 0.3)

	plt.setp(rnaFrac_axis.xaxis.get_majorticklabels(), rotation = 20, fontsize = FONT_SIZE)
	plt.setp(proteinFrac_axis.xaxis.get_majorticklabels(), rotation = 20, fontsize = FONT_SIZE)
	plt.setp(dryMassInit_axis.xaxis.get_majorticklabels(), rotation = 20, fontsize = FONT_SIZE)
	plt.setp(dryMassFinal_axis.xaxis.get_majorticklabels(), rotation = 20, fontsize = FONT_SIZE)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


def getMaxTime(allCells):
	maxTime = 0
	for simDir in allCells:
		simOutDir = os.path.join(simDir, "simOut")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		maxTime = np.max([maxTime, time.size])
	return maxTime

def getSingleValue(allCells, tableName, colName, maxTime):
	allCellsData = np.ones((allCells.size, maxTime), np.float64) * np.nan

	for idx, simDir in enumerate(allCells):
		simOutDir = os.path.join(simDir, "simOut")
		value = TableReader(os.path.join(simOutDir, tableName)).readColumn(colName)
		allCellsData[idx, :value.size] = value

	return allCellsData

def removeNanReshape(a):
	return a[np.logical_not(np.isnan(a))]

def getExpectedComposition(doubling_time):
	# return np.ones(6), 100. * units.fg

	from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
	raw_data = KnowledgeBaseEcoli()

	from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
	sim_data = fitSimData_1(raw_data, doubling_time = doubling_time)

	subMasses = sim_data.mass.avgCellSubMass
	proteinMass = subMasses['proteinMass'].asNumber(units.fg) / sim_data.mass.avgCellToInitialCellConvFactor
	rnaMass = subMasses['rnaMass'].asNumber(units.fg) / sim_data.mass.avgCellToInitialCellConvFactor

	initialMass = sim_data.mass.avgCellDryMassInit.asNumber(units.fg)

	return proteinMass, rnaMass, initialMass

def downSample(vector, factor):
	return vector[:, np.arange(start = 0, stop = vector.shape[1], step = factor)]

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
