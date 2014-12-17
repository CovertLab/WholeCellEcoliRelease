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
	
	# Load parameters
	measured_synthetase_rates = unfitKb.trna_synthetase_rates
	fit_synthetase_rates = mostfitKb.trna_synthetase_rates.asNumber()
	amino_acid_labels = mostfitKb.aa_trna_groups.keys()

	synthetase_mean = mostfitKb.synthetase_counts
	synthetase_variance = mostfitKb.synthetase_variance
	initial_aa_poly_rate = mostfitKb.initial_aa_polymerization_rate.asNumber()
	min_turnover_rates = mostfitKb.minimum_trna_synthetase_rates.asNumber()

	# Plotting
	plt.figure(figsize = (8.5, 11))

	ax = plt.subplot(3,1,1)

	plt.plot(fit_synthetase_rates, '_', label='fit', markeredgewidth=3)
	plt.plot(measured_synthetase_rates, 'x', label='Jakubowski et al. 1984', markeredgewidth=2)
	plt.plot(min_turnover_rates, '_', label='min fit', markeredgewidth=3)

	plt.xlabel('amino acid')
	plt.ylabel('1/s')
	plt.title('trna synthetase turnover rates')
	plt.xticks(np.arange(0,21), amino_acid_labels)

	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.8])

	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, prop={'size':9})

	plt.subplot(3,1,2)
	plt.bar(np.arange(synthetase_mean.size), synthetase_mean, yerr=synthetase_variance)
	plt.xlabel('amino acid')
	plt.ylabel('mean initial synthetase count')
	plt.xticks(np.arange(0,21) + 0.4, amino_acid_labels)

	plt.subplot(3,1,3)
	plt.bar(np.arange(initial_aa_poly_rate.size), initial_aa_poly_rate)
	plt.xlabel('amino acid')
	plt.ylabel('initial rate of aa polymerization (aa/s)')
	plt.xticks(np.arange(0,21) + 0.4, amino_acid_labels)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


if __name__ == "__main__":
	defaultKBDirectory = wholecell.utils.constants.SERIALIZED_KB_DIR

	parser = argparse.ArgumentParser()
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbDirectory", help = "KB file name", type = str, default = defaultKBDirectory)

	args = parser.parse_args().__dict__

	main(args["plotOutDir"], args["plotOutFileName"], args["kbDirectory"])