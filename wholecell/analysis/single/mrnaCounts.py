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

import tables
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.utils.knowledgebase_fixture_manager import loadKnowledgeBase

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = loadKnowledgeBase(simOutDir)

	isMRna = kb.rnaData["isMRna"]

	rnaIds = kb.rnaData["id"][isMRna]

	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:

		names = bulkMoleculesFile.root.names
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules

		moleculeIds = names.moleculeIDs.read()

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]

	# avgCounts = rnaCountsBulk.mean(0)

	# relativeCounts = avgCounts / avgCounts.sum()

	# relativeCounts = rnaCountsBulk[-1, :] / rnaCountsBulk[-1, :].sum()

	plt.figure(figsize = (8.5, 11))

	counts = rnaCountsBulk[-1, :]

	expectedCountsArbitrary = kb.rnaExpression['expression'][isMRna].magnitude

	expectedCounts = expectedCountsArbitrary/expectedCountsArbitrary.sum() * counts.sum()

	maxLine = 1.1 * max(expectedCounts.max(), counts.max())
	plt.plot([0, maxLine], [0, maxLine], '--r')	
	plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

	plt.xlabel("Expected RNA count (scaled to total)")
	plt.ylabel("Actual RNA count (at final time step)")

	# plt.show()

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])
