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
from matplotlib_venn import venn3, venn3_circles

data = np.array([1181, 382, 1832, 80, 877, 0, 1]) #order is 100, 010, 110, 001, 101, 011, 111

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fig = plt.figure()
	v = venn3(subsets = data, set_labels = ["0 < Freq. < 1", "Freq. = 1", "Freq. = 0"])
	v.get_patch_by_id("100").set_color("green")
	v.get_patch_by_id("010").set_color("blue")
	v.get_patch_by_id("001").set_color("red")
	v.get_patch_by_id("110").set_color("cyan")
	v.get_patch_by_id("101").set_color("yellow")
	v.get_patch_by_id("111").set_color("white")

	ids = ["100", "010", "001", "110", "101", "111"]

	for i in ids:
		v.get_patch_by_id(i).set_edgecolor("none")
		v.get_patch_by_id(i).set_alpha(0.8)

		if i == "100":
			v.get_label_by_id(i).set_text("1,181")

		if i == "110":
			v.get_label_by_id(i).set_text("1,832")

	setLabels = v.set_labels
	setLabels[0].set_color("green")
	setLabels[1].set_color("blue")
	setLabels[2].set_color("red")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	for i in ids:
		v.get_label_by_id(i).set_text("")
	for i in setLabels:
		i.set_text("")

	exportFigure(plt, plotOutDir, plotOutFileName + "_clean", metadata)
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