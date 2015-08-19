#!/usr/bin/env python
"""
Plots time-step effects on polymerization

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	nAvogadro = kb.constants.nAvogadro
	ribosomeSubunitIds = []
	ribosomeSubunitIds.extend(kb.moleculeGroups.s50_fullComplex)
	ribosomeSubunitIds.extend(kb.moleculeGroups.s30_fullComplex)
	ribosomeSubunitIds.extend(kb.moleculeGroups.s50_proteinComplexes)
	ribosomeSubunitIds.extend(kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitIds'])
	ribosomeSubunitIds.extend(kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitIds'])
	ribosomeSubunitMasses = kb.getter.getMass(ribosomeSubunitIds)
	mass70s = (kb.getter.getMass(kb.moleculeGroups.s50_fullComplex) + kb.getter.getMass(kb.moleculeGroups.s30_fullComplex))[0]

	elongationRate = float(kb.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s))

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load ribosome data
	ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
	ribosomesTerminated = ribosomeDataFile.readColumn("didTerminate")
	ribosomesInitialized = ribosomeDataFile.readColumn("didInitialize")
	ribosomeTerminationLoss = ribosomeDataFile.readColumn("terminationLoss")
	ribosomeDataFile.close()

	# Load RNAP data
	rnapDataFile = TableReader(os.path.join(simOutDir, "RnapData"))
	rnapsTerminated = rnapDataFile.readColumn("didTerminate")
	rnapsInitialized = rnapDataFile.readColumn("didInitialize")
	rnapTerminationLoss = rnapDataFile.readColumn("terminationLoss")
	rnapDataFile.close()

	# Load count data for s30 proteins, rRNA, and final 30S complex
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	# Get indexes
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)

	# Load data
	ribosomeSubunitCounts = bulkMolecules.readColumn("counts")[:, ribosomeSubunitIndexes]

	bulkMolecules.close()

	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

	uniqueMoleculeCounts.close()

	# Calculate statistics
	totalRibosome = (activeRibosome + np.min(ribosomeSubunitCounts))
	totalRibosomeCapacity = totalRibosome * elongationRate

	freeSubunitMass = (ribosomeSubunitMasses * ribosomeSubunitCounts / nAvogadro).asNumber(units.fg).sum(axis = 1)
	activeRibosomeMass = (mass70s * activeRibosome / nAvogadro).asNumber(units.fg)
	totalRibosomeMass = freeSubunitMass + activeRibosomeMass
	massFractionActive = activeRibosomeMass / totalRibosomeMass

	plt.figure(figsize = (8.5, 11))

	ribosomeInit_axis = plt.subplot(8,1,1)
	ribosomeInit_axis.plot(time / 60., ribosomesInitialized, label="Number of ribosomes initialized", linewidth=2, color='b')
	ribosomeInit_axis.set_ylabel("ribosomes")
	ribosomeInit_axis.legend()#ncol=2)

	ribosomeTerm_axis = plt.subplot(8,1,2)
	ribosomeTerm_axis.plot(time / 60., ribosomesTerminated, label="Number of ribosomes terminated", linewidth=2, color='r')
	ribosomeTerm_axis.set_ylabel("ribosomes")
	ribosomeTerm_axis.legend()#ncol=2)

	ribosomeTermLoss_axis = plt.subplot(8,1,3)
	ribosomeTermLoss_axis.plot(time / 60., ribosomeTerminationLoss, label="Lost capacity due to termination", linewidth=2, color='r')
	ribosomeTermLoss_axis.set_ylabel("lost aa capacity")
	ribosomeTermLoss_axis.legend()#ncol=2)

	rnapInit_axis = plt.subplot(8,1,4)
	rnapInit_axis.plot(time / 60., rnapsInitialized, label="Number of rnap initialized", linewidth=2, color='r')
	rnapInit_axis.set_ylabel("rnap")
	rnapInit_axis.legend()#ncol=2)

	rnapTerm_axis = plt.subplot(8,1,5)
	rnapTerm_axis.plot(time / 60., rnapsTerminated, label="Number of rnap terminated", linewidth=2, color='r')
	rnapTerm_axis.set_ylabel("rnap")
	rnapTerm_axis.legend()#ncol=2)

	ribosomeTermLoss_axis = plt.subplot(8,1,6)
	ribosomeTermLoss_axis.plot(time / 60., rnapTerminationLoss, label="Lost capacity due to termination", linewidth=2, color='r')
	ribosomeTermLoss_axis.set_ylabel("lost ntp capacity")
	ribosomeTermLoss_axis.legend()#ncol=2)

	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
