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

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	ribosomeSubunitIds = [kb.moleculeGroups.s50_fullComplex[0], kb.moleculeGroups.s30_fullComplex[0]]
	elongationRate = float(kb.constants.ribosomeElongationRate.asNumber(units.aa / units.s))

	# Load ribosome data
	massFile = TableReader(os.path.join(simOutDir, "RibosomeData"))

	actualElongations = massFile.readColumn("actualElongations")
	expectedElongations_recorded = massFile.readColumn("expectedElongations")
	time = massFile.readColumn("time")

	massFile.close()

	# Load count data for s30 proteins, rRNA, and final 30S complex
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	# Get indexes
	moleculeIds = bulkMolecules.readAttribute("moleculeIDs")
	ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)

	# Load data
	ribosomeSubunitCounts = bulkMolecules.readColumn("counts")[:, ribosomeSubunitIndexes]

	bulkMolecules.close()

	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

	uniqueMoleculeCounts.close()

	# Calculate total ribosome elongation capacity
	activeRibosomeCapacity = activeRibosome * elongationRate
	totalRibosomeCapacity = (activeRibosome + np.min(ribosomeSubunitCounts)) * elongationRate

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	ribosomeCapacity_axis = plt.subplot(2,1,1)
	ribosomeCapacity_axis.plot(time / 60., totalRibosomeCapacity, label="Total ribosome capacity", linewidth=2, color='b')
	ribosomeCapacity_axis.plot(time / 60., actualElongations, label="Actual elongations", linewidth=2, color='r')
	ribosomeCapacity_axis.set_ylabel("amino acids polymerized")
	ribosomeCapacity_axis.legend(ncol=2)



	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

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
