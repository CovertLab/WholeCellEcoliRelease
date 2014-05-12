#!/usr/bin/env python
"""
Plot ribosome counts

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

	RIBOSOME_RNA_IDS = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]
	ribosomeRnaIndexes = np.array([moleculeIds.index(rRnaId) for rRnaId in RIBOSOME_RNA_IDS], np.int)
	ribosomeRnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, ribosomeRnaIndexes]

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "UniqueMolecules.hdf"))

	ribosome = h.root.activeRibosome
	isActive = ribosome.col("_entryState") == 1
	time = ribosome.col("_time")
	nActive = np.bincount(time[isActive])

	time = np.unique(time)

	plt.figure(figsize = (8.5, 11))

	plt.plot(time / 60, nActive[time])
	plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[time[0]], 2 * nActive[time[0]]], "r--")
	plt.xlabel("Time (min)")
	plt.ylabel("Counts")
	plt.title("Active Ribosomes Final:Initial = %0.2f" % (nActive[time[-1]] / float(nActive[time[0]])))

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

	h.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])