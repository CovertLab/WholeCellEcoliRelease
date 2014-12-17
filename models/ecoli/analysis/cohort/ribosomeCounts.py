#!/usr/bin/env python

"""
Plot ribosome counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/19/2014
"""

from __future__ import division

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.containers.unique_molecules_data import UniqueMoleculesData

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	simOutSubDirs = sorted([
		os.path.join(simOutDir, item)
		for item in os.listdir(simOutDir)
		if os.path.isdir(os.path.join(simOutDir, item)) and item != "kb"
		])

	plt.figure(figsize = (8.5, 11))

	for simOutSubDir in simOutSubDirs:
		print simOutSubDir

		h = tables.open_file(os.path.join(simOutSubDir, "BulkMolecules.hdf"))

		names = h.root.names
		bulkMolecules = h.root.BulkMolecules

		moleculeIds = names.moleculeIDs.read()

		RIBOSOME_RNA_IDS = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]
		ribosomeRnaIndexes = np.array([moleculeIds.index(rRnaId) for rRnaId in RIBOSOME_RNA_IDS], np.int)
		ribosomeRnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, ribosomeRnaIndexes]

		h.close()

		h = tables.open_file(os.path.join(simOutSubDir, "UniqueMoleculeCounts.hdf"))

		uniqueMoleculeCounts = h.root.UniqueMoleculeCounts
		ribosomeIndex = uniqueMoleculeCounts.attrs.uniqueMoleculeIds.index("activeRibosome")
		time = uniqueMoleculeCounts.col("time")
		nActive = uniqueMoleculeCounts.col("uniqueMoleculeCounts")[:, ribosomeIndex]

		h.close()

		plt.plot(time / 60, nActive)
		plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("Counts")
	
	plt.title("Active Ribosomes")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

	h.close()

if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])