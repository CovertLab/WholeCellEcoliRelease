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

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir)

	# Get all cells
	allDir = ap.getAll()

	#plt.figure(figsize = (8.5, 11))

	compoundsPlotted = set()
	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")
		#initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		counts = TableReader(os.path.join(simOutDir, "BulkMolecules")).readColumn("counts")
		allNames = TableReader(os.path.join(simOutDir, "BulkMolecules")).readAttribute('objectNames')
				
		compoundNames = []
		nonZeroCounts = counts.T[np.any(counts.T, axis = 1)]
		for idx, compound in enumerate(nonZeroCounts):
			if idx > 100:
				break
			compartment = allNames[idx][-3:]
			if compartment == '[c]':
				compoundNames.append(allNames[idx][:-3])
				compoundsPlotted.add(allNames[idx])
				plt.plot(time / 60. / 60., compound)

		plt.legend(compoundNames)
			
		# # set axes to size that shows all generations
		# cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
		# if cellCycleTime > currentMaxTime:
		# 	currentMaxTime = cellCycleTime

		# axesList[idx].set_xlim(0, currentMaxTime*int(metadata["total_gens"])*1.1)
		# axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")


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
