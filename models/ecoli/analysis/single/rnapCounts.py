#!/usr/bin/env python
"""
Plot RNA polymerase counts and counts of mRNA precursors

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnapId = "APORNAP-CPLX[c]"
	rnapIndex = moleculeIds.index(rnapId)
	rnapCountsBulk = bulkMolecules.readColumn("counts")[:, rnapIndex]

	RNAP_RNA_IDS = ["EG10893_RNA[c]", "EG10894_RNA[c]", "EG10895_RNA[c]", "EG10896_RNA[c]"]
	rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNAP_RNA_IDS], np.int)
	rnapRnaCounts = bulkMolecules.readColumn("counts")[:, rnapRnaIndexes]

	bulkMolecules.close()

	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

	rnapIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	nActive = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnapIndex]

	uniqueMoleculeCounts.close()

	plt.figure(figsize = (8.5, 11))

	plt.subplot(5, 1, 1)

	plt.plot(time / 60., nActive + rnapCountsBulk)
	plt.xlabel("Time (min)")
	plt.ylabel("Protein Counts")
	plt.title("RNA Polymerase")

	for subplotIdx in xrange(2, 6):
		rnapRnaCountsIdx = subplotIdx - 2

		plt.subplot(5, 1, subplotIdx)

		plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA counts")
		plt.title(RNAP_RNA_IDS[rnapRnaCountsIdx])

	plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


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
