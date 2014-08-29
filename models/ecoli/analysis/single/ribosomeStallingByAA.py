#!/usr/bin/env python
"""
Plot ribosome stalling on a per amino acid basis

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/27/2014
"""

from __future__ import division

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	with tables.open_file(os.path.join(simOutDir, "RibosomeStalling.hdf")) as h5file:
		timeStep = h5file.root.RibosomeStalling.col("timeStep")
		aaCountInSequence = h5file.root.RibosomeStalling.col("aaCountInSequence")
		aaCounts = h5file.root.RibosomeStalling.col("aaCounts")
		trnaCapacity = h5file.root.RibosomeStalling.col("trnasCapacity")
		synthetaseCapacity = h5file.root.RibosomeStalling.col("synthetaseCapacity")

	aaLimitation = -1 * (aaCountInSequence - aaCounts).clip(min = 0)
	trnaCapacityLimitation = -1 * (aaCountInSequence - trnaCapacity).clip(min = 0)
	synthetaseCapacityLimitation = -1 * (aaCountInSequence - synthetaseCapacity).clip(min = 0)

	aaExcess = -1 * (aaCountInSequence - aaCounts).clip(max = 0)
	trnaCapacityExcess = -1 * (aaCountInSequence - trnaCapacity).clip(max = 0)
	synthetaseCapacityExcess = -1 * (aaCountInSequence - synthetaseCapacity).clip(max = 0)

	kb = cPickle.load(open(kbFile, "rb"))
	amino_acid_labels = kb.aa_trna_groups.keys()

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(21):
		plt.subplot(6, 4, idx + 1)

		plt.plot(timeStep / 60., aaLimitation[:,idx])
		plt.plot(timeStep / 60., trnaCapacityLimitation[:,idx])
		plt.plot(timeStep / 60., synthetaseCapacityLimitation[:,idx])
		plt.plot(timeStep / 60., aaExcess[:,idx])
		plt.plot(timeStep / 60., trnaCapacityExcess[:,idx])
		plt.plot(timeStep / 60., synthetaseCapacityExcess[:,idx])

		plt.title(amino_acid_labels[idx])

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)
	plt.tight_layout()
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
