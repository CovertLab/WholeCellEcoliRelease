#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
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

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "EvaluationTime.hdf"))

	evalTimes = h.root.EvaluationTime

	stateNames = evalTimes.attrs.stateNames
	processNames = evalTimes.attrs.processNames

	time = evalTimes.col("time")

	updateQueries_times = evalTimes.col("updateQueries_times")
	partition_times = evalTimes.col("partition_times")
	merge_times = evalTimes.col("merge_times")
	calculateRequest_times = evalTimes.col("calculateRequest_times")
	evolveState_times = evalTimes.col("evolveState_times")
	# updateQueries_total = evalTimes.col("updateQueries_total")
	# partition_total = evalTimes.col("partition_total")
	# merge_total = evalTimes.col("merge_total")
	# calculateRequest_total = evalTimes.col("calculateRequest_total")
	# evolveState_total = evalTimes.col("evolveState_total")

	h.close()

	plt.figure(figsize = (8.5, 11))

	plt.subplot(3, 2, 1)

	plt.semilogy(time / 60, updateQueries_times * 1000)
	plt.grid(True, which = "major")
	plt.xlabel("Simulation time (min)")
	plt.ylabel("Evaluation time (ms)")
	plt.title("State.updateQueries")
	plt.legend(stateNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

	plt.subplot(3, 2, 2)

	plt.semilogy(time / 60, calculateRequest_times * 1000)
	plt.grid(True, which = "major")
	plt.xlabel("Simulation time (min)")
	plt.ylabel("Evaluation time (ms)")
	plt.title("Process.calculateRequest")
	plt.legend(processNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

	plt.subplot(3, 2, 3)

	plt.semilogy(time / 60, partition_times * 1000)
	plt.grid(True, which = "major")
	plt.xlabel("Simulation time (min)")
	plt.ylabel("Evaluation time (ms)")
	plt.title("State.partition")
	plt.legend(stateNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

	plt.subplot(3, 2, 4)

	plt.semilogy(time / 60, evolveState_times * 1000)
	plt.grid(True, which = "major")
	plt.xlabel("Simulation time (min)")
	plt.ylabel("Evaluation time (ms)")
	plt.title("Process.evolveState")
	plt.legend(processNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

	plt.subplot(3, 2, 5)

	plt.semilogy(time / 60, merge_times * 1000)
	plt.grid(True, which = "major")
	plt.xlabel("Simulation time (min)")
	plt.ylabel("Evaluation time (ms)")
	plt.title("State.merge")
	plt.legend(stateNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

	plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

	h.close()

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
