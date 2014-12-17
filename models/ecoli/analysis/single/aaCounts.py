#!/usr/bin/env python
"""
Plot amino acid counts

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

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf"))

	names = h.root.names

	moleculeIds = names.moleculeIDs.read()

	AA_IDS = [
		"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]",
		"CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]",
		"HIS-L[c]", "ILE-L[c]", "LEU-L[c]", "LYS-L[c]",
		"MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]",
		"THR-L[c]", "TRP-L[c]", "TYR-L[c]", "SEC-L[c]",
		"VAL-L[c]"
		]
	aaIndexes = np.array([moleculeIds.index(aaId) for aaId in AA_IDS], np.int)
	bulkMolecules = h.root.BulkMolecules
	aaCounts = bulkMolecules.read(0, None, 1, "counts")[:, aaIndexes]

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))
	table = h.root.Mass
	time = np.array([x["time"] for x in table.iterrows()])
	h.close()

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(21):

		plt.subplot(6, 4, idx + 1)

		plt.plot(time / 60., aaCounts[:, idx], linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("Counts")
		plt.title(AA_IDS[idx])

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
