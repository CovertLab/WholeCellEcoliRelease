#!/usr/bin/env python
"""
Plot amino acid counts

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/13/2014
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

from wholecell.containers.unique_molecules_data import bundleByFieldValue

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	TABLE_FILE_REPLACE = TableReader(os.path.join(simOutDir, "ReplicationForkPosition"))
	data = h.root.ReplicationForkPosition.read()
	TABLE_FILE_REPLACE.close()

	plt.figure(figsize = (8.5, 11))

	legendEntries = []

	for dnaPolyId, dnaPolyHistory in bundleByFieldValue(data, "_uniqueId", ["_timeStep", "chromosomeLocation"]):
		timeMinutes = dnaPolyHistory["_timeStep"] / 60
		position = dnaPolyHistory["chromosomeLocation"]

		plt.plot(timeMinutes, position, linewidth = 2)

		legendEntries.append(dnaPolyId)

	plt.xlabel("Time (min)")
	plt.ylabel("Replication fork position (nt)")

	plt.legend(legendEntries, loc = "best")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

if __name__ == "__main__":
	# defaultKBFile = os.path.join(
	# 		wholecell.utils.constants.SERIALIZED_KB_DIR,
	# 		wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
	# 		)

	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	# parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])

	print "Still need to fix this, nontrivial - John"
