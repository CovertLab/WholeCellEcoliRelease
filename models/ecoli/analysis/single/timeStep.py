#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	timeSteps = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	meanTimeStep = np.mean(timeSteps)
	maxTimeStep = np.max(timeSteps)
	minTimeStep = np.min(timeSteps)


	plt.figure(figsize = (8.5, 11))

	plt.subplot(2,1,1)
	plt.title("Simulation Time Steps")
	plt.plot(timeSteps)
	plt.xlabel("Increment")
	plt.ylabel("timeStep (s)")

	plt.subplot(2,1,2)
	plt.figtext(.15, .035, "Mean = %.2f, Max = %.2f, Min = %.2f" % (meanTimeStep, maxTimeStep, minTimeStep))
	plt.plot(time,timeSteps)
	plt.xlabel("Cell time (s)")
	plt.ylabel("timeStep (s)")


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

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
