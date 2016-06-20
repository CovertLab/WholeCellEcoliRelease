#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS

BURN_IN_SECONDS = 500
DISABLED = True

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	if DISABLED:
		print "Currently disabled because it requires too much memory."
		return

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# Get all cells
	allDir = ap.get_cells()

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		counts = TableReader(os.path.join(simOutDir, "BulkMolecules")).readColumn("counts")
		countsToMolar = TableReader(os.path.join(simOutDir, "EnzymeKinetics")).readColumn("countsToMolar")
		allNames = TableReader(os.path.join(simOutDir, "BulkMolecules")).readAttribute('objectNames')

		compoundNames = []
		nonZeroCounts = counts.T[np.any(counts.T, axis = 1)]
		for idx, counts in enumerate(nonZeroCounts):
			if (counts[BURN_IN_SECONDS:] > 0).sum() > 100:
				compartment = allNames[idx][-3:]
				compoundNames.append(allNames[idx][:20])
				concentrations = (counts * countsToMolar)
				if time[0] < 1:
					concentrations[:BURN_IN_SECONDS] = np.mean(concentrations[BURN_IN_SECONDS:])
				plt.plot(time / 60., concentrations / np.mean(concentrations))

		# plt.legend(compoundNames, fontsize=5)
		plt.title("Protein Concentrations")
		plt.xlabel("Time (min)")
		plt.ylabel("Mean-normalized concentration")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "fixme")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
