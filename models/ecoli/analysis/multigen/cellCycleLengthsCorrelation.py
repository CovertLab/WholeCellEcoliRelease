#!/usr/bin/env python
"""
Plots correlation of cell cycle lengths between each generation and the next.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/24/2015
"""
import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(seedOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir)


	# Get all cells
	allDir = ap.getAll()

	cellCycleLengths = []
	generations = []
	for idx, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		cellCycleLengths.append((time[-1] - time[0]) / 60. / 60.)
		generations.append(idx)

	differences = [1, 2, 3]
	for difference in differences:
		xArray = cellCycleLengths[0:-difference]
		yArray = cellCycleLengths[difference:]
		coefficient, pValue = scipy.stats.pearsonr(xArray, yArray)
		
		try:
			plt.subplot(1, len(differences), difference, sharey=axis)
		except NameError:
			axis = plt.subplot(1, len(differences), difference)

		plt.scatter(xArray, yArray)
		plt.xlabel('Cell cycle length of Generation n')
		plt.ylabel('Cell cycle length of Generation n+%i' % (difference))
		plt.text(np.amin(xArray), np.amax(yArray),'r = %.4f' % (coefficient))
		plt.xticks([np.amin(xArray), (np.amin(xArray) + np.amax(xArray)) / 2. ,np.amax(xArray)])
		plt.yticks([np.amin(yArray), (np.amin(yArray) + np.amax(yArray)) / 2. ,np.amax(yArray)])
		
		
	plt.title('Correlation of cell cycle lengths')	
	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

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
