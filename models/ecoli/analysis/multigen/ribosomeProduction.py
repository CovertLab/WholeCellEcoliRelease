#!/usr/bin/env python

from __future__ import division

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

	## Get expected doubling time ##
	expected_doubling_time = sim_data.doubling_time

	fig = plt.figure()
	fig.set_size_inches(10,12)

	for gen, simDir in enumerate(firstCellLineage):
		simOutDir = os.path.join(simDir, "simOut")

		## Mass growth rate ##
		time, growthRate = getMassData(simDir, ["instantaniousGrowthRate"])
		timeStep = units.s * TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
		time = units.s * time
		growthRate = (1 / units.s) * growthRate
		doublingTime = 1 / growthRate * np.log(2)

		## Calculate ribosomal rna doubling times ##
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

		rrn16S_doubling_time = 1 / ( (1 / timeStep) * (rrn16S_produced / rrn16s_count) ) * np.log(2)
		rrn23S_doubling_time = 1 / ( (1 / timeStep) * (rrn23S_produced / rrn23s_count) ) * np.log(2)
		rrn5S_doubling_time = 1 / ( (1 / timeStep) * (rrn5S_produced / rrn5s_count) ) * np.log(2)
		
		rrn16S_doubling_time[rrn16S_doubling_time.asNumber() == np.inf] = np.nan * units.s
		rrn23S_doubling_time[rrn23S_doubling_time.asNumber() == np.inf] = np.nan * units.s
		rrn5S_doubling_time[rrn5S_doubling_time.asNumber() == np.inf] = np.nan * units.s

		## Calculate ribosomal rna actual initiation probabilities ##
		rrn16S_init_prob = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn16S_init_prob")
		rrn23S_init_prob = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn23S_init_prob")
		rrn5S_init_prob = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn5S_init_prob")

		idx_16s = np.zeros(len(sim_data.moleculeGroups.s30_16sRRNA), dtype=np.int)
		for idx, id16s in enumerate(sim_data.moleculeGroups.s30_16sRRNA):
			idx_16s[idx] = np.where(sim_data.process.transcription.rnaData['id'] == id16s)[0][0]

		idx_23s = np.zeros(len(sim_data.moleculeGroups.s50_23sRRNA), dtype=np.int)
		for idx, id23s in enumerate(sim_data.moleculeGroups.s50_23sRRNA):
			idx_23s[idx] = np.where(sim_data.process.transcription.rnaData['id'] == id23s)[0][0]

		idx_5s = np.zeros(len(sim_data.moleculeGroups.s50_5sRRNA), dtype=np.int)
		for idx, id5s in enumerate(sim_data.moleculeGroups.s50_5sRRNA):
			idx_5s[idx] = np.where(sim_data.process.transcription.rnaData['id'] == id5s)[0][0]

		rrn16s_fit_init_prob = sim_data.process.transcription.rnaData['synthProb'][idx_16s].sum()
		rrn23s_fit_init_prob = sim_data.process.transcription.rnaData['synthProb'][idx_23s].sum()
		rrn5s_fit_init_prob = sim_data.process.transcription.rnaData['synthProb'][idx_5s].sum()

		## Calculated expected multinomial variance ##
		total_rna_init = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("total_rna_init")
		rrn16s_exp_variance = total_rna_init * rrn16s_fit_init_prob * (1 - rrn16s_fit_init_prob)
		rrn23s_exp_variance = total_rna_init * rrn23s_fit_init_prob * (1 - rrn23s_fit_init_prob)
		rrn5s_exp_variance = total_rna_init * rrn5s_fit_init_prob * (1 - rrn5s_fit_init_prob)

		rrn16s_actual_variance = rrn16s_fit_init_prob.std() ** 2
		rrn23s_actual_variance = rrn23s_fit_init_prob.std() ** 2
		rrn5s_actual_variance = rrn5s_fit_init_prob.std() ** 2

		## Plotting ##
		gs = gridspec.GridSpec(7, 3)

		# Plot growth rate
		ax1 = plt.subplot(gs[0,:2])
		avgDoublingTime = doublingTime[1:].asNumber(units.min).mean()
		stdDoublingTime = doublingTime[1:].asNumber(units.min).std()
		ax1.plot(time.asNumber(units.min), doublingTime.asNumber(units.min))
		ax1.plot(time.asNumber(units.min), expected_doubling_time.asNumber(units.min) * np.ones(time.asNumber().size), linestyle='--')
		ax1.set_ylim([avgDoublingTime - 2*stdDoublingTime, avgDoublingTime + 2*stdDoublingTime])
		ax1.set_ylabel("Doubling\ntime (min)")
		ax1.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')

		ax1_1 = plt.subplot(gs[0,2])
		hist_doublingTime = removeNanReshape(doublingTime.asNumber(units.min))
		nbins = np.ceil(np.sqrt(hist_doublingTime.size))
		ax1_1.hist(hist_doublingTime, nbins)

		##

		ax2 = plt.subplot(gs[1,:2])
		ax2.plot(time.asNumber(units.min), rrn16S_doubling_time.asNumber(units.min))
		ax2.plot(time.asNumber(units.min), expected_doubling_time.asNumber(units.min) * np.ones(time.asNumber().size), linestyle='--')
		ax2.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax2.set_ylabel("rrn 16S\ndoubling time")

		ax2_1 = plt.subplot(gs[1,2])
		hist_rrn16S_doubling_time = removeNanReshape(rrn16S_doubling_time.asNumber(units.min))
		nbins = np.ceil(np.sqrt(hist_rrn16S_doubling_time.size))
		ax2_1.hist(hist_rrn16S_doubling_time, nbins)

		##

		ax3 = plt.subplot(gs[2,:2])
		ax3.plot(time.asNumber(units.min), rrn23S_doubling_time.asNumber(units.min))
		ax3.plot(time.asNumber(units.min), expected_doubling_time.asNumber(units.min) * np.ones(time.asNumber().size), linestyle='--')
		ax3.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax3.set_ylabel("rrn 23S\ndoubling time")

		ax3_1 = plt.subplot(gs[2,2])
		hist_rrn23S_doubling_time = removeNanReshape(rrn23S_doubling_time.asNumber(units.min))
		nbins = np.ceil(np.sqrt(hist_rrn23S_doubling_time.size))
		ax3_1.hist(hist_rrn23S_doubling_time, nbins)

		##

		ax4 = plt.subplot(gs[3,:2])
		ax4.plot(time.asNumber(units.min), rrn5S_doubling_time.asNumber(units.min))
		ax4.plot(time.asNumber(units.min), expected_doubling_time.asNumber(units.min) * np.ones(time.asNumber().size), linestyle='--')
		ax4.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax4.set_ylabel("rrn 5S\ndoubling time")

		ax4_1 = plt.subplot(gs[3,2])
		hist_rrn5S_doubling_time = removeNanReshape(rrn5S_doubling_time.asNumber(units.min))
		nbins = np.ceil(np.sqrt(hist_rrn5S_doubling_time.size))
		ax4_1.hist(hist_rrn5S_doubling_time, nbins)

		##

		ax5 = plt.subplot(gs[4,:2])
		ax5.plot(time.asNumber(units.min), rrn16S_init_prob)
		ax5.plot(time.asNumber(units.min), rrn16s_fit_init_prob * np.ones(time.asNumber().size), linestyle='--')
		ax5.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax5.set_ylabel("rrn 16S\ninit prob")

		ax5_1 = plt.subplot(gs[4,2])
		hist_rrn16S_init_prob = removeNanReshape(rrn16S_init_prob / total_rna_init)
		nbins = np.ceil(np.sqrt(hist_rrn16S_init_prob.size))
		ax5_1.hist(hist_rrn16S_init_prob, nbins)

		##

		ax6 = plt.subplot(gs[5,:2])
		ax6.plot(time.asNumber(units.min), rrn23S_init_prob)
		ax6.plot(time.asNumber(units.min), rrn23s_fit_init_prob * np.ones(time.asNumber().size), linestyle='--')
		ax6.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax6.set_ylabel("rrn 23S\ninit prob")

		ax6_1 = plt.subplot(gs[5,2])
		hist_rrn23S_init_prob = removeNanReshape(rrn23S_init_prob / total_rna_init)
		nbins = np.ceil(np.sqrt(hist_rrn23S_init_prob.size))
		ax6_1.hist(hist_rrn23S_init_prob, nbins)

		##

		ax7 = plt.subplot(gs[6,:2])
		ax7.plot(time.asNumber(units.min), rrn5S_init_prob)
		ax7.plot(time.asNumber(units.min), rrn5s_fit_init_prob * np.ones(time.asNumber().size), linestyle='--')
		ax7.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax7.set_ylabel("rrn 5S\ninit prob")

		ax7_1 = plt.subplot(gs[6,2])
		hist_rrn5S_init_prob = removeNanReshape(rrn5S_init_prob / total_rna_init)
		nbins = np.ceil(np.sqrt(hist_rrn5S_init_prob.size))
		ax7_1.hist(hist_rrn5S_init_prob, nbins)

		ax7.set_xlabel("Time (min)")

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

def removeNanReshape(a):
	return a[np.logical_not(np.isnan(a))]

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