#!/usr/bin/env python
"""
Plot ribosome stalling

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/22/2014
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
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

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ribosomeData = TableReader(os.path.join(simOutDir, "RibosomeData"))

	timeStep = ribosomeData.readColumn("timeStep")
	# stallingRateTotal = ribosomeData.readColumn("stallingRateTotal")
	# stallingRateMean = ribosomeData.readColumn("stallingRateMean")
	# stallingRateStd = ribosomeData.readColumn("stallingRateStd")
	fractionStalled = ribosomeData.readColumn("fractionStalled")
	aaCountInSequence = ribosomeData.readColumn("aaCountInSequence")
	aaCounts = ribosomeData.readColumn("aaCounts")
	trnaCapacity = ribosomeData.readColumn("trnasCapacity")
	synthetaseCapacity = ribosomeData.readColumn("synthetaseCapacity")

	ribosomeData.close()

	aaLimitation = -1 * (aaCountInSequence - aaCounts).clip(min = 0).sum(axis = 1)
	trnaCapacityLimitation = -1 * (aaCountInSequence - trnaCapacity).clip(min = 0).sum(axis = 1)
	synthetaseCapacityLimitation = -1 * (aaCountInSequence - synthetaseCapacity).clip(min = 0).sum(axis = 1)

	aaExcess = -1 * (aaCountInSequence - aaCounts).clip(max = 0).sum(axis = 1)
	trnaCapacityExcess = -1 * (aaCountInSequence - trnaCapacity).clip(max = 0).sum(axis = 1)
	synthetaseCapacityExcess = -1 * (aaCountInSequence - synthetaseCapacity).clip(max = 0).sum(axis = 1)

	plt.figure(figsize = (8.5, 11))
	plt.subplot(2,1,1)
	plt.plot(timeStep / 60, fractionStalled)

	plt.xlabel("Time (min)")
	plt.ylabel("Fraction of ribosomes stalled")

	plt.subplot(2,1,2)

	plt.plot(timeStep / 60, aaLimitation, '--', label = 'aa limit')
	plt.plot(timeStep / 60, trnaCapacityLimitation, '--', label = 'trna limit')
	plt.plot(timeStep / 60, synthetaseCapacityLimitation, '--', label = 'synthetase limit')

	plt.plot(timeStep / 60, aaExcess, label = 'aa excess')
	plt.plot(timeStep / 60, trnaCapacityExcess, label = 'trna excess')
	plt.plot(timeStep / 60, synthetaseCapacityExcess, label = 'synthetase excess')
	plt.legend(prop={'size':7})
	plt.xlabel("Time (min)")
	plt.ylabel("Magnitude of capacity/demand mismatch (elongations)")

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


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
