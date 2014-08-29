#!/usr/bin/env python

"""
justAnalyzeFitKb.py

Creates plots base don unfit and most-fit KB's for fitting analysis
"""

import cPickle
import os
import argparse

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants

def main(plotOutDir, plotOutFileName, kbDirectory):
	unfitKbFileName = os.join(kbDirectory, wholecell.utils.constants.SERIALIZED_KB_UNFIT_FILENAME)
	mostfitKbFileName = os.join(kbDirectory, wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME)

	unfitKb = cPickle.load(open(unfitKbFileName, "rb"))
	mostfitKb = cPickle.load(open(mostfitKbFileName, "rb"))

	plt.figure(figsize = (8.5, 11))

	measured_synthetase_rates = unfitKb.

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


if __name__ == "__main__":
	defaultKBDirectory = wholecell.utils.constants.SERIALIZED_KB_DIR

	parser = argparse.ArgumentParser()
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbDirectory", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["plotOutDir"], args["plotOutFileName"], args["kbDirectory"])