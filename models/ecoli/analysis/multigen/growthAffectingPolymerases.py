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

		## RNAP counts and statistics ##

		# Get free counts
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnapId = "APORNAP-CPLX[c]"
		rnapIndex = moleculeIds.index(rnapId)
		rnapCountsBulk = bulkMolecules.readColumn("counts")[:, rnapIndex]
		bulkMolecules.close()

		# Get active counts
		uniqueMolecules = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		rnapIndex = uniqueMolecules.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		rnapCountsActive = uniqueMolecules.readColumn("uniqueMoleculeCounts")[:, rnapIndex]
		uniqueMolecules.close()

		# Calculate statistics
		totalRnap = rnapCountsBulk + rnapCountsActive
		fractionRnapActive = rnapCountsActive / (rnapCountsActive + rnapCountsBulk)

		## Ribosome counts and statistics ##

		# Get ids for 30S and 50S subunits
		proteinIds30S = sim_data.moleculeGroups.s30_proteins
		rnaIds30S = [sim_data.process.translation.monomerData['rnaId'][np.where(sim_data.process.translation.monomerData['id'] == pid)[0][0]] for pid in proteinIds30S]
		rRnaIds30S = sim_data.moleculeGroups.s30_16sRRNA
		complexIds30S = [sim_data.moleculeGroups.s30_fullComplex[0]]

		proteinIds50S = sim_data.moleculeGroups.s50_proteins
		rnaIds50S = [sim_data.process.translation.monomerData['rnaId'][np.where(sim_data.process.translation.monomerData['id'] == pid)[0][0]] for pid in proteinIds50S]
		rRnaIds50S = sim_data.moleculeGroups.s50_23sRRNA
		rRnaIds50S.extend(sim_data.moleculeGroups.s50_5sRRNA)
		complexIds50S = [sim_data.moleculeGroups.s50_fullComplex[0]]

		# Get indexes for 30S and 50S subunits based on ids
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		proteinIndexes30S = np.array([moleculeIds.index(protein) for protein in proteinIds30S], np.int)
		rnaIndexes30S = np.array([moleculeIds.index(rna) for rna in rnaIds30S], np.int)
		rRnaIndexes30S = np.array([moleculeIds.index(rRna) for rRna in rRnaIds30S], np.int)
		complexIndexes30S = np.array([moleculeIds.index(comp) for comp in complexIds30S], np.int)

		proteinIndexes50S = np.array([moleculeIds.index(protein) for protein in proteinIds50S], np.int)
		rnaIndexes50S = np.array([moleculeIds.index(rna) for rna in rnaIds50S], np.int)
		rRnaIndexes50S = np.array([moleculeIds.index(rRna) for rRna in rRnaIds50S], np.int)
		complexIndexes50S = np.array([moleculeIds.index(comp) for comp in complexIds50S], np.int)

		# Get counts of 30S and 50S mRNA, rProteins, rRNA, and full complex counts
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		freeProteinCounts30S = bulkMolecules.readColumn("counts")[:, proteinIndexes30S]
		rnaCounts30S = bulkMolecules.readColumn("counts")[:, rnaIndexes30S]
		freeRRnaCounts30S = bulkMolecules.readColumn("counts")[:, rRnaIndexes30S]
		complexCounts30S = bulkMolecules.readColumn("counts")[:, complexIndexes30S]
		
		freeProteinCounts50S = bulkMolecules.readColumn("counts")[:, proteinIndexes50S]
		rnaCounts50S = bulkMolecules.readColumn("counts")[:, rnaIndexes50S]
		freeRRnaCounts50S = bulkMolecules.readColumn("counts")[:, rRnaIndexes50S]
		complexCounts50S = bulkMolecules.readColumn("counts")[:, complexIndexes50S]

		bulkMolecules.close()

		# Get active ribosome counts
		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		uniqueMoleculeCounts.close()

		# Get elongation rate data
		ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
		actualElongations = ribosomeDataFile.readColumn("actualElongations")
		ribosomeDataFile.close()

		# Calculate statistics
		rProteinCounts = np.hstack((freeProteinCounts50S, freeProteinCounts30S))
		limitingCountThreshold = 15
		limitingRProteinIdxs = np.unique(np.where(rProteinCounts < limitingCountThreshold)[1])
		limitingRProteinCounts = rProteinCounts[:, limitingRProteinIdxs] # Get traces of limiting rProteins

		rRnaCounts = np.hstack((freeRRnaCounts50S, freeRRnaCounts30S))
		rRnaCounts = rRnaCounts[:, np.unique(np.where(rRnaCounts > 0)[1])] # Get only non-zero for all time counts

		counts30S = complexCounts30S
		counts50S = complexCounts50S

		ribosomeCounts = activeRibosome

		effectiveElongationRate = actualElongations / ribosomeCounts
		extraRibosomes = (ribosomeCounts - actualElongations / 21.) / (actualElongations / 21.) * 100

		## Calculate statistics involving ribosomes and RNAP ##
		ratioRNAPtoRibosome = totalRnap.astype(np.float) / ribosomeCounts.astype(np.float)

		## Plotting ##

		gs = gridspec.GridSpec(9, 3)

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

		# Plot RNAP active fraction
		ax2 = plt.subplot(gs[1,:2])
		ax2.plot(time.asNumber(units.min), fractionRnapActive)
		ax2.plot(time.asNumber(units.min), sim_data.growthRateParameters.fractionActiveRnap * np.ones(time.asNumber().size), linestyle='--')
		ax2.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax2.set_ylabel("Fraction active\nRNAP")

		ax2_1 = plt.subplot(gs[1,2])
		hist_fractionRnapActive = removeNanReshape(fractionRnapActive)
		nbins = np.ceil(np.sqrt(hist_fractionRnapActive.size))
		ax2_1.hist(hist_fractionRnapActive, nbins)

		# Plot RNAP active and total counts
		ax3 = plt.subplot(gs[2,:2])
		ax3.plot(time.asNumber(units.min), totalRnap)
		ax3.plot(time.asNumber(units.min), rnapCountsActive)
		ax3.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax3.set_ylabel("Total & active\nRNAP")

		# Plot limiting rProtein counts
		ax4 = plt.subplot(gs[3,:2])
		if limitingRProteinCounts.size > 0:
			ax4.plot(time.asNumber(units.min), limitingRProteinCounts)
		ax4.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax4.set_ylim([0, 100])
		ax4.set_ylabel("Limiting rProtein\ncounts")

		# Plot rRNA counts
		ax5 = plt.subplot(gs[4,:2])
		ax5.plot(time.asNumber(units.min), rRnaCounts)
		ax5.set_ylim([0, 200])
		ax5.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax5.set_ylabel("rRNA\ncounts")

		# Plot 30S and 50S counts
		ax6 = plt.subplot(gs[5,:2])
		ax6.plot(time.asNumber(units.min), counts30S)
		ax6.plot(time.asNumber(units.min), counts50S)
		ax6.set_ylim([0, np.max([counts30S[2:].max(), counts50S[2:].max()])])
		ax6.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax6.set_ylabel("30S & 50S\ncounts")

		# Plot ribosome counts
		ax7 = plt.subplot(gs[6,:2])
		ax7.plot(time.asNumber(units.min), ribosomeCounts)
		ax7.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax7.set_ylabel("Ribosome counts")

		# Plot RNAP:ribosome ratio
		ax8 = plt.subplot(gs[7,:2])
		ax8.plot(time.asNumber(units.min), ratioRNAPtoRibosome)
		ax8.plot(time.asNumber(units.min), ratioRNAPtoRibosome.mean() * np.ones(time.asNumber().size), linestyle='--')
		ax8.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax8.set_ylabel("RNAP:Ribosome\ncounts")

		# Plot number of "extra" ribosomes
		ax9 = plt.subplot(gs[8,:2])
		ax9.plot(time.asNumber(units.min), extraRibosomes)
		ax9.axvline(x = time.asNumber(units.min).max(), linewidth=2, color='k', linestyle='--')
		ax9.set_ylim([0, 100])
		ax9.set_ylabel("% extra\nribosomes")

		ax9.set_xlabel("Time (min)")

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
