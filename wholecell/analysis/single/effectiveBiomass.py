#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import argparse
import os
import cPickle

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as sch

import wholecell.utils.constants

FLUX_UNITS = "M/s"

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	with tables.open_file(os.path.join(simOutDir, "FBAResults.hdf")) as h5file:
		timeStep = h5file.root.FBAResults.col("timeStep")
		outputFluxes = h5file.root.FBAResults.col("outputFluxes")

		names = h5file.root.names
		outputMoleculeIDs = np.array(names.outputMoleculeIDs.read())

	fig = plt.figure(figsize = (30, 15))

	normalized = (outputFluxes / np.max(np.abs(outputFluxes[20:, :]), 0)).transpose()

	linkage = sch.linkage(normalized, method = "centroid")
	dendro = sch.dendrogram(linkage, orientation="right", no_plot=True)
	index = dendro["leaves"]

	ax = fig.gca()

	ax.matshow(
		np.fmax(np.fmin((normalized[index, :]+1)/2, +1), 0), # colormap expects values from 0 to 1
		cmap = plt.get_cmap("coolwarm")
		)

	ax.set_yticks(np.arange(len(index)))
	ax.set_yticklabels(outputMoleculeIDs[np.array(index)], size = 5)

	plt.title("Relative FBA production rates (red = production, blue = consumption)")

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


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
