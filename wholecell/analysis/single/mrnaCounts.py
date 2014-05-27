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
# matplotlib.use("Agg")
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

	relativeCounts = rnaCountsBulk[-1, :] / rnaCountsBulk[-1, :].sum()

	expectedCounts = kb.rnaExpression['expression'][isMRna].magnitude

	expectedRelativeCounts = expectedCounts/expectedCounts.sum()

	plt.plot(expectedRelativeCounts, relativeCounts, '.')

	plt.show()

	import ipdb; ipdb.set_trace()

	# plt.figure(figsize = (8.5, 11))

	# plt.plot(time / 60, nActive)
	# plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")
	# plt.xlabel("Time (min)")
	# plt.ylabel("Counts")
	# plt.title("Active Ribosomes Final:Initial = %0.2f" % (nActive[-1] / float(nActive[0])))

	# plt.savefig(os.path.join(plotOutDir, plotOutFileName))


if __name__ == "__main__":
	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])

	# TODO: reinstate this analysis
	pass
