#!/usr/bin/env python
"""
Plot mRNA counts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	isMRna = kb.process.transcription.rnaData["isMRna"]

	rnaIds = kb.process.transcription.rnaData["id"][isMRna]

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]

	# avgCounts = rnaCountsBulk.mean(0)

	# relativeCounts = avgCounts / avgCounts.sum()

	# relativeCounts = rnaCountsBulk[-1, :] / rnaCountsBulk[-1, :].sum()

	bulkMolecules.close()

	plt.figure(figsize = (8.5, 11))

	counts = rnaCountsBulk[-1, :]

	expectedCountsArbitrary = kb.process.transcription.rnaData['expression'][isMRna]

	expectedCounts = expectedCountsArbitrary/expectedCountsArbitrary.sum() * counts.sum()

	maxLine = 1.1 * max(expectedCounts.max(), counts.max())
	plt.plot([0, maxLine], [0, maxLine], '--r')
	plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

	plt.xlabel("Expected RNA count (scaled to total)")
	plt.ylabel("Actual RNA count (at final time step)")

	# plt.show()

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
