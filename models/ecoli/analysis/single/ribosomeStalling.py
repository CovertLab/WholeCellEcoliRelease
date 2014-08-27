#!/usr/bin/env python
"""
Plot ribosome stalling

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/22/2014
"""

from __future__ import division

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

	with tables.open_file(os.path.join(simOutDir, "RibosomeStalling.hdf")) as h5file:
		timeStep = h5file.root.RibosomeStalling.col("timeStep")
		# stallingRateTotal = h5file.root.RibosomeStalling.col("stallingRateTotal")
		# stallingRateMean = h5file.root.RibosomeStalling.col("stallingRateMean")
		# stallingRateStd = h5file.root.RibosomeStalling.col("stallingRateStd")
		fractionStalled = h5file.root.RibosomeStalling.col("fractionStalled")
		aaCountInSequence = h5file.root.RibosomeStalling.col("aaCountInSequence")
		aaCounts = h5file.root.RibosomeStalling.col("aaCounts")
		trnasCapacity = h5file.root.RibosomeStalling.col("trnasCapacity")
		synthetaseCapacity = h5file.root.RibosomeStalling.col("synthetaseCapacity")


	import ipdb; ipdb.set_trace()
	
	plt.figure(figsize = (8.5, 11))
	plt.subplot(2,1,1)
	plt.plot(timeStep / 60, fractionStalled)

	plt.xlabel("Time (min)")
	plt.ylabel("Fraction of ribosomes stalled")


	plt.subplot(2,1,2)

	plt.xlabel("Time (min)")
	plt.ylabel("Stalling caused by fraction")

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
