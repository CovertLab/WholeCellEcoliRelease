"""
Plots correlation of cell cycle lengths between each generation and the next.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/24/2015
"""

from __future__ import absolute_import

import os

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)


		# Get all cells
		allDir = ap.get_cells()

		plt.figure(figsize = (10, 10))

		cellCycleLengths = []
		for idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

			cellCycleLengths.append((time[-1] - time[0]) / 60. / 60.)

		totalTimeAverage = np.average(cellCycleLengths)
		totalGens = int(metadata["total_gens"])

		if totalGens > 8:
			differences = [1,2]
			for value in np.arange(4, totalGens, totalGens//4):
				differences.append(value)
		elif totalGens < 4:
			return
		else:
			differences = [1,2]

		for idx, difference in enumerate(differences):
			xArray = cellCycleLengths[0:-difference]
			yArray = cellCycleLengths[difference:]
			coefficient, pValue = scipy.stats.pearsonr(xArray, yArray)

			xMin = np.amin(xArray)
			xMax = np.amax(xArray)
			yMin = np.amin(yArray)
			yMax = np.amax(yArray)
			yMax = np.fmax(yMax, 1.1 * yMin)

			try:
				plt.subplot(1, len(differences), idx+1, sharey=axis)
			except NameError:
				axis = plt.subplot(1, len(differences), idx+1)

			plt.scatter(xArray, yArray)
			plt.ylabel('Cell cycle length of Generation n+%i' % (difference,))
			plt.text(xMin, yMax,'r = %.4f' % (coefficient,))
			plt.xticks([xMin, (xMin + xMax) / 2., xMax], 
				[str(round(xMin,4)), str(round((xMin + xMax) / 2.,4)), str(round(xMax,4))],
				rotation='vertical')
			plt.yticks([yMin, (yMin + yMax) / 2., yMax])


		plt.xlabel('Cell cycle length of Generation n')
		plt.suptitle('Correlation of cell cycle lengths')
		plt.subplots_adjust(hspace = 0.9, wspace = 0.8, bottom = 0.2)


		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
