#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from wholecell.utils.sparkline import whitePadSparklineAxis

SHUFFLE_VARIANT_TAG = "ShuffleParams"
PLACE_HOLDER = -1

FONT_SIZE=9
trim = 0.05


def mm2inch(value):
	return value * 0.0393701

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):

	if metadata is not None and SHUFFLE_VARIANT_TAG not in metadata["variant"]:
		print "This plot only runs for variants where parameters are shuffled."
		return

	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)
	control_sim = ap.get_cells(variant = [0])
	variant_cells = ap.get_cells(variant = range(1, ap.n_variant))

	doublingTimes = []
	for simDir in variant_cells:
		try:
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")

			doublingTimes.append( (time.max() - initialTime) / 60. )
		except:
			continue

	doublingTimes = np.array(doublingTimes)

	controlDoublingTime = None
	for simDir in control_sim:

		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")

		controlDoublingTime = (time.max() - initialTime) / 60.



	fig = plt.figure()
	fig.set_figwidth(5)
	fig.set_figheight(5)
	ax = plt.subplot(1, 1, 1)
	ax.hist(doublingTimes, np.sqrt(doublingTimes.size))
	ax.axvline(controlDoublingTime, color = "k", linestyle = "dashed", linewidth = 2)

	ax.set_xlabel("Cell Division Time (min)")
	ax.set_title("Mean: %0.3g     Std: %0.3g     Control: %0.3g" % (doublingTimes.mean(), doublingTimes.std(), controlDoublingTime))

	axes_list = [ax]

	for a in axes_list:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in a.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

	whitePadSparklineAxis(ax)

	plt.subplots_adjust(bottom = 0.2, wspace=0.3)

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
