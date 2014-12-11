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

	with tables.open_file(os.path.join(simOutDir, "RibosomeData.hdf")) as h5file:
		timeStep = h5file.root.RibosomeData.col("timeStep")
		aaCountInSequence = h5file.root.RibosomeData.col("aaCountInSequence")
		aaCounts = h5file.root.RibosomeData.col("aaCounts")
		trnaCapacity = h5file.root.RibosomeData.col("trnasCapacity")
		synthetaseCapacity = h5file.root.RibosomeData.col("synthetaseCapacity")

	aaCapacity = -1 * (aaCountInSequence - aaCounts)
	trnaCapacity = -1 * (aaCountInSequence - trnaCapacity)
	synthetaseCapacity = -1 * (aaCountInSequence - synthetaseCapacity)

	kb = cPickle.load(open(kbFile, "rb"))

	amino_acid_labels = [
		"ALA", "ARG", "ASN", "ASP",
		"CYS", "GLU", "GLN", "GLY",
		"HIS", "ILE", "LEU", "LYS",
		"MET", "PHE", "PRO", "SER",
		"THR", "TRP", "TYR", "SEC",
		"VAL"
		]

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(21):
		plt.subplot(6, 4, idx + 1)

		plt.plot(timeStep / 60., aaCapacity[:,idx], linewidth = 2,label = 'aa limit')
		plt.plot(timeStep / 60., trnaCapacity[:,idx], '--', label = 'trna limit')
		plt.plot(timeStep / 60., synthetaseCapacity[:,idx], label = 'synthetase limit')

		plt.title(amino_acid_labels[idx])
	plt.legend(bbox_to_anchor=(1.02, 0.5, 4., .102), loc=5, ncol=2, mode="expand")

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
