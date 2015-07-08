#!/usr/bin/env python
"""
Plot protein monomer counts

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
from wholecell.utils import units

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of proteins from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	proteinIds = kb.process.translation.monomerData["id"]

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

	proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]

	bulkMolecules.close()

	# avgCounts = proteinCountsBulk.mean(0)

	# relativeCounts = avgCounts / avgCounts.sum()

	counts = proteinCountsBulk[-1, :]

	expectedCountsArbitrary = (
		kb.process.transcription.rnaData["expression"][kb.relation.rnaIndexToMonomerMapping] /
		(np.log(2) / kb.doubling_time.asNumber(units.s) + kb.process.translation.monomerData["degRate"].asNumber(1/units.s))
		) * counts.sum()

	expectedCountsRelative = expectedCountsArbitrary / expectedCountsArbitrary.sum()

	expectedCounts = expectedCountsRelative * counts.sum()

	plt.figure(figsize = (8.5, 11))

	maxLine = 1.1 * max(expectedCounts.max(), counts.max())
	plt.plot([0, maxLine], [0, maxLine], '--r')
	plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

	plt.xlabel("Expected protein count (scaled to total)")
	plt.ylabel("Actual protein count (at final time step)")

	# plt.show()

	# plt.plot(time / 60, nActive)
	# plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")
	# plt.xlabel("Time (min)")
	# plt.ylabel("Counts")
	# plt.title("Active Ribosomes Final:Initial = %0.2f" % (nActive[-1] / float(nActive[0])))

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
