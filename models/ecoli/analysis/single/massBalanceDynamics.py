#!/usr/bin/env python

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

THRESHOLD = 1e-5

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	mass = TableReader(os.path.join(simOutDir, "Mass"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = (
		TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		) / 60.

	relProcessMassDifferences = mass.readColumn("relProcessMassDifferences")

	relProcessMassDifferences[relProcessMassDifferences == 0] = np.nan

	processNames = mass.readAttribute("processNames")

	mass.close()

	plt.figure(figsize = (8.5, 11))

	n_processes = len(processNames)

	n_cols = int(np.sqrt(n_processes))
	n_rows = int(np.ceil(n_processes/n_cols))

	axis = [time.min(), time.max(), np.nanmin(np.abs(relProcessMassDifferences)), np.nanmax(np.abs(relProcessMassDifferences))]

	for i, processName in enumerate(processNames):
		plt.subplot(n_rows, n_cols, i+1)

		series = relProcessMassDifferences[:, i].copy()
		t = time.copy()

		t = t[series != 0]
		series = series[series != 0]

		if np.any(series > 0):
			colors = np.zeros((series.size, 3), np.float64)

			colors[np.abs(series) > THRESHOLD, :] = (1.0, 0.0, 0.0)

			# markers = ["^" if value > 0 else "v" for value in series]

			plt.scatter(
				t,
				np.abs(series),
				color = colors,
				marker = '.'
				# marker = markers
				)

		plt.hlines(THRESHOLD,  axis[0], axis[1], "r", "dashed")

		plt.title(processName)

		# plt.xlabel("time (min)")
		# plt.ylabel("mass diff ($|(m_f - m_i)/(m_i)|$)")

		plt.yscale("log")

		plt.axis(axis)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
