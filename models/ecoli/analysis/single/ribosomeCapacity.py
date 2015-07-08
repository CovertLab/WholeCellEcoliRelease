#!/usr/bin/env python
"""
Plots ribosome capacity

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/20/2014
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
from wholecell.utils.sparkline import sparklineAxis, setAxisMaxMinY

FONT = {
		'size'	:	8
		}

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

	elongationRate = float(kb.constants.ribosomeElongationRate.asNumber(units.aa / units.s))

	# Load ribosome data
	ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))

	actualElongations = ribosomeDataFile.readColumn("actualElongations")
	expectedElongations_recorded = ribosomeDataFile.readColumn("expectedElongations")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	ribosomeDataFile.close()

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
	totalRibosome = (activeRibosome + ribosomeSubunitCounts.min(axis=1))
	totalRibosomeCapacity = totalRibosome * elongationRate

	freeSubunitMass = (ribosomeSubunitMasses * ribosomeSubunitCounts / nAvogadro).asNumber(units.fg).sum(axis = 1)
	activeRibosomeMass = (mass70s * activeRibosome / nAvogadro).asNumber(units.fg)
	totalRibosomeMass = freeSubunitMass + activeRibosomeMass
	massFractionActive = activeRibosomeMass / totalRibosomeMass

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	ribosomeCapacity_axis = plt.subplot(6,1,1)
	ribosomeCapacity_axis.plot(time / 60., totalRibosomeCapacity, label="Theoretical total ribosome capacity", linewidth=2, color='b')
	ribosomeCapacity_axis.plot(time / 60., actualElongations, label="Actual elongations", linewidth=2, color='r')
	ribosomeCapacity_axis.set_ylabel("Amino acids polymerized")
	ribosomeCapacity_axis.legend(ncol=2)

	activeRibosomeCapacity_axis = plt.subplot(6,1,2)
	activeRibosomeCapacity_axis.plot(time / 60., activeRibosome * elongationRate, label="Theoretical active ribosome capacity", linewidth=2, color='b')
	activeRibosomeCapacity_axis.plot(time / 60., actualElongations, label="Actual elongations", linewidth=2, color='r')
	activeRibosomeCapacity_axis.set_ylabel("Amino acids polymerized")
	activeRibosomeCapacity_axis.legend(ncol=2)

	inactiveRibosomeCapacity_axis = plt.subplot(6,1,3)
	inactiveRibosomeCapacity_axis.plot(time / 60., ribosomeSubunitCounts.min(axis=1) * elongationRate, label="Theoretical inactive ribosome capacity", linewidth=2, color='b')
	inactiveRibosomeCapacity_axis.set_ylabel("Amino acids polymerized")
	inactiveRibosomeCapacity_axis.legend(ncol=2)

	fractionalCapacity_axis = plt.subplot(6,1,4)
	fractionalCapacity_axis.plot(time / 60., actualElongations / totalRibosomeCapacity, label="Fraction of ribosome capacity used", linewidth=2, color='k')
	fractionalCapacity_axis.set_ylabel("Fraction of ribosome capacity used")
	fractionalCapacity_axis.set_yticks(np.arange(0., 1.05, 0.05))
	#fractionalCapacity_axis.get_yaxis().grid(b=True, which='major', color='b', linestyle='--')
	fractionalCapacity_axis.grid(b=True, which='major', color='b', linestyle='--')

	effectiveElongationRate_axis = plt.subplot(6,1,5)
	effectiveElongationRate_axis.plot(time / 60., actualElongations / activeRibosome, label="Effective elongation rate", linewidth=2, color='k')
	effectiveElongationRate_axis.set_ylabel("Effective elongation rate (aa/s/ribosome)")

	fractionActive_axis = plt.subplot(6,1,6)
	fractionActive_axis.plot(time / 60., massFractionActive, label="Mass fraction active", linewidth=2, color='k')
	fractionActive_axis.set_ylabel("Mass fraction of active ribosomes")
	fractionActive_axis.set_yticks(np.arange(0., 1.1, 0.1))

	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.6)

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
