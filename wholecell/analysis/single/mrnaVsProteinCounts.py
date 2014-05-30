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

	rnaIds = kb.rnaData["id"][kb.rnaIndexToMonomerMapping]

	proteinIds = kb.monomerData["id"]

	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:

		names = bulkMoleculesFile.root.names
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules

		moleculeIds = names.moleculeIDs.read()

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]

		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

		proteinCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, proteinIndexes]

	relativeMRnaCounts = rnaCountsBulk[-1, :] #/ rnaCountsBulk[-1, :].sum()
	relativeProteinCounts = proteinCountsBulk[-1, :] #/ proteinCountsBulk[-1, :].sum()

	plt.figure(figsize = (8.5, 11))

	plt.plot(relativeMRnaCounts, relativeProteinCounts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

	plt.xlabel("RNA count (at final time step)")
	plt.ylabel("Protein count (at final time step)")

	# plt.show()

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])
