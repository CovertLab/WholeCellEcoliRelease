#!/usr/bin/env python
"""
Plot mass fractions

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

def main(simOutDir, plotOutDir, plotOutFileName):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf"))

	names = h.root.names
	bulkMolecules = h.root.BulkMolecules

	moleculeIds = names.moleculeIDs.read()
	RNAP_IDS = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]
	rnapIndexes = np.array([moleculeIds.index(rnapId) for rnapId in RNAP_IDS], np.int)
	rnapCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnapIndexes]

	RNAP_RNA_IDS = ["EG10893_RNA[c]", "EG10894_RNA[c]", "EG10895_RNA[c]", "EG10896_RNA[c]"]
	rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNAP_RNA_IDS], np.int)
	rnapRnaCounts = bulkMolecules.read(0, None, 1, "counts")[:, rnapRnaIndexes]

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "UniqueMolecules.hdf"))

	rnaPol = h.root.activeRnaPoly
	isActive = rnaPol.col("_entryState") == 1
	time = rnaPol.col("_time")
	nActive = np.bincount(time[isActive])

	time = np.unique(time)

	plt.figure(figsize = (8.5, 11))

	plt.subplot(5, 1, 1)

	plt.plot(time / 60., nActive[time] + np.min(rnapCountsBulk, axis = 1)[1:])
	plt.xlabel("Time (min)")
	plt.ylabel("Protein Counts")
	plt.title("RNA Polymerase")

	for subplotIdx in xrange(2, 6):
		rnapRnaCountsIdx = subplotIdx - 2
	
		plt.subplot(5, 1, subplotIdx)

		plt.plot(time / 60., rnapRnaCounts[1:, rnapRnaCountsIdx])
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA counts")
		plt.title(RNAP_RNA_IDS[rnapRnaCountsIdx])

	plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

	h.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])