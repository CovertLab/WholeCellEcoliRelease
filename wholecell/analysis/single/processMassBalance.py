#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2014
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

THRESHOLD = 1e-9 # roughly, one hundredth of a hydrogen atom

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	with tables.open_file(os.path.join(simOutDir, "Mass.hdf")) as h5file:
		table = h5file.root.Mass

		time = table.col("time")
		processMassDifferences = table.col("processMassDifferences")

		processNames = table.attrs.processNames

	avgProcessMassDifferences = np.abs(processMassDifferences).sum(0) / len(time)

	index = np.arange(len(processNames))
	width = 1

	plt.figure(figsize = (8.5, 11))

	axes = plt.axes()

	r1 = axes.barh(index, avgProcessMassDifferences * (avgProcessMassDifferences > THRESHOLD), width, log = True, color = (0.9, 0.2, 0.2))
	r2 = axes.barh(index, avgProcessMassDifferences * (avgProcessMassDifferences <= THRESHOLD), width, log = True, color = (0.2, 0.2, 0.9))

	axes.set_yticks(index+width/2)
	axes.set_yticklabels(processNames) #, rotation = -45)

	axes.plot([THRESHOLD, THRESHOLD], [index[0], index[-1]+width], 'k--', linewidth=3)

	# axes.legend((r1, r2), ("Positive", "Negative"))

	plt.xlabel("Mass difference (fg)")

	plt.title("Average absolute change in mass by individual processes")

	plt.tight_layout()
	plt.grid(True, which = "major")

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
