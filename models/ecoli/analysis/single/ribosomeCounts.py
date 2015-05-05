#!/usr/bin/env python
"""
Plot ribosome counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

# from wholecell.containers.unique_molecules_data import UniqueMoleculesData

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	RIBOSOME_RNA_IDS = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]
	ribosomeRnaIndexes = np.array([moleculeIds.index(rRnaId) for rRnaId in RIBOSOME_RNA_IDS], np.int)
	ribosomeRnaCountsBulk = bulkMolecules.readColumn("counts")[:, ribosomeRnaIndexes]

	bulkMolecules.close()

	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	nActive = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

	uniqueMoleculeCounts.close()

	plt.figure(figsize = (8.5, 11))

	plt.plot(time / 60, nActive)
	plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")
	plt.xlabel("Time (min)")
	plt.ylabel("Counts")
	plt.title("Active Ribosomes Final:Initial = %0.2f" % (nActive[-1] / float(nActive[0])))

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
