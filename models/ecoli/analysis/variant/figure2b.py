#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from wholecell.utils.sparkline import whitePadSparklineAxis

FONT_SIZE=9

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):

	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	if ap.n_generation == 1:
		print "Need more data to create plot"
		return

	fig = plt.figure()
	fig.set_figwidth(15)
	fig.set_figheight(5)

	doublingTimeVariants = [44, 100, 22]

	for varIdx in range(ap.n_variant):

		if varIdx == 0:
			plotIdx = 1
			gen = [2,3]
		elif varIdx == 1:
			plotIdx = 0
			gen = [2,3]
		elif varIdx == 2:
			plotIdx = 2
			gen = [6,7]
		else:
			continue

		initial_masses = np.zeros(0)
		final_masses = np.zeros(0)

		all_cells = ap.get_cells(generation=[2,3], variant=[varIdx])
		if len(all_cells) == 0:
			continue

		doublingTimes = np.zeros(len(all_cells))
		for idx, simDir in enumerate(all_cells):
			try:
				simOutDir = os.path.join(simDir, "simOut")
				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			except:
				print 'Error reading data for %s' % (simDir)

			doublingTimes[idx] = (time[-1] - time[0]) / 60.

		bins = 16
		ax = plt.subplot2grid((1, 3), (0, plotIdx))
		ax.hist(doublingTimes, bins)
		ax.axvline(x = doublingTimeVariants[varIdx], color = "r", linestyle = "--")

		ax.set_title("%i min" % (doublingTimeVariants[varIdx]), fontsize = FONT_SIZE)

		ax.set_xlabel("Doubling Time (min)", fontsize = FONT_SIZE)

		plt.subplots_adjust(bottom = 0.2)

		whitePadSparklineAxis(ax)

		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

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
