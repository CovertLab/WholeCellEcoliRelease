#!/usr/bin/env python
"""
Plots ribosome capacity

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/20/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

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
	ribosomeSubunitIds = [kb.s50_fullComplex, kb.s30_fullComplex]
	elongationRate = float(kb.ribosomeElongationRate.asNumber(units.aa / units.s))
	
	# Load ribosome data
	with tables.open_file(os.path.join(simOutDir, "RibosomeData.hdf")) as massFile:
		table = massFile.root.RibosomeData
		actualElongations = table.col("actualElongations")
		expectedElongations_recorded = table.col("expectedElongations")
		time = table.col("time")

	# Load count data for s30 proteins, rRNA, and final 30S complex
	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules

		# Get indexes
		moleculeIds = bulkMoleculesFile.root.names.moleculeIDs.read()
		ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)

		# Load data
		ribosomeSubunitCounts = bulkMolecules.read(0, None, 1, "counts")[:, ribosomeSubunitIndexes]

	with tables.open_file(os.path.join(simOutDir, "UniqueMoleculeCounts.hdf")) as uniqueMoleculesFile:
		uniqueMoleculeCounts = uniqueMoleculesFile.root.UniqueMoleculeCounts
		ribosomeIndex = uniqueMoleculeCounts.attrs.uniqueMoleculeIds.index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.col("uniqueMoleculeCounts")[:, ribosomeIndex]

	# Calculate total ribosome elongation capacity
	activeRibosomeCapacity = activeRibosome * elongationRate
	totalRibosomeCapacity = (activeRibosome + np.min(ribosomeSubunitCounts)) * elongationRate

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	totalRibosomeCapacity_axis = plt.subplot(2,1,1)
	sparklineAxis(totalRibosomeCapacity_axis, time / 60., totalRibosomeCapacity, 'left', '-', 'b')
	setAxisMaxMinY(totalRibosomeCapacity_axis, totalRibosomeCapacity)
	totalRibosomeCapacity_axis.set_ylabel("Total ribosome capacity (aa)")

	actualElongations_axis = totalRibosomeCapacity_axis.twinx()
	sparklineAxis(actualElongations_axis, time / 60., actualElongations, 'right', '-', 'r')
	setAxisMaxMinY(actualElongations_axis, totalRibosomeCapacity)
	actualElongations_axis.set_ylabel("Actual elongations (aa)")

	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

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
