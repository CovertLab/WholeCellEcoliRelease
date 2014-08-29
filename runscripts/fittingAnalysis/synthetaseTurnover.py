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
import numpy as np

import wholecell.utils.constants

def main(plotOutDir, plotOutFileName, kbDirectory):
	unfitKbFileName = os.path.join(kbDirectory, wholecell.utils.constants.SERIALIZED_KB_UNFIT_FILENAME)
	mostfitKbFileName = os.path.join(kbDirectory, wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME)

	unfitKb = cPickle.load(open(unfitKbFileName, "rb"))
	mostfitKb = cPickle.load(open(mostfitKbFileName, "rb"))

	measured_synthetase_rates = unfitKb.trna_synthetase_rates
	fit_synthetase_rates = mostfitKb.trna_synthetase_rates
	amino_acid_labels = mostfitKb.aa_trna_groups.keys()


	# Plotting
	plt.figure(figsize = (8.5, 11))

	plt.plot(measured_synthetase_rates, 'x', label='Jakubowski et al. 1984')
	plt.plot(fit_synthetase_rates, 'o', label='fit')

	plt.xlabel('amino acid')
	plt.ylabel('1/s')
	plt.title('trna synthetase turnover rates')
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=5, ncol=2, mode="expand")
	plt.xticks(np.arange(0,21), amino_acid_labels)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))
	#plt.show()


if __name__ == "__main__":
	defaultKBDirectory = wholecell.utils.constants.SERIALIZED_KB_DIR

	parser = argparse.ArgumentParser()
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbDirectory", help = "KB file name", type = str, default = defaultKBDirectory)

	args = parser.parse_args().__dict__

	main(args["plotOutDir"], args["plotOutFileName"], args["kbDirectory"])