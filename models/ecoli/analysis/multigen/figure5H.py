#!/usr/bin/env python
"""
Plots Figure 5D. Data is obtained from an analysis of compendium of environment shifts.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2017
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import wholecell.utils.constants
from matplotlib_venn import venn3

data = [1143, 330, 1920, 79, 854, 2, 25] #order is 100, 010, 110, 001, 101, 011, 111
dataArea = [1143, 330, 1920, 79, 854, 500, 25]

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fig = plt.figure()
	v = venn3(subsets = dataArea, set_labels = ["0 < Freq. < 1", "Freq. = 1", "Freq. = 0"])
	ids = ["100", "010", "001", "110", "101", "011", "111"]

	for i, color in zip(ids, ["b", "r", "y", "purple", "green", "orange", "white"]):
		v.get_patch_by_id(i).set_color(color)

	for label, color in zip(v.set_labels, ["b", "r", "y"]):
		label.set_color(color)

	v.get_label_by_id("011").set_text("%s" % data[-2])

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	for i in ids:
		v.get_label_by_id(i).set_text("")

	for label in v.set_labels:
		label.set_text("")

	plt.savefig(os.path.join(plotOutDir, plotOutFileName + "__clean.pdf"))
	plt.close()

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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])