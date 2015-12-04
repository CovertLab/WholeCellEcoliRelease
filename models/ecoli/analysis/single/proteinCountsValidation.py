#!/usr/bin/env python
"""
Compare protein counts to Wisniewski 2014 data set

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/3/2015
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
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)


	sim_data = cPickle.load(open(simDataFile, "rb"))
	validation_data = cPickle.load(open(validationDataFile, "rb"))

	ids_complex = sim_data.process.complexation.moleculeNames
	ids_eqBind = sim_data.process.equilibrium.moleculeNames
	ids_monomer = sim_data.process.translation.monomerData["id"].tolist()
	ids_protein = sorted(set(ids_complex + ids_eqBind + ids_monomer))
	bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
	view_complex = bulkContainer.countsView(ids_complex)
	view_eqBind = bulkContainer.countsView(ids_eqBind)
	view_monomer = bulkContainer.countsView(ids_monomer)
	view_validation = bulkContainer.countsView(validation_data.protein.wisniewski2014Data["monomerId"].tolist())


	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")

	proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], np.int)
	proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]
	bulkMolecules.close()

	bulkContainer.countsIs(proteinCountsBulk.mean(axis = 0))

	# TODO: IMPORTANT
	# TODO: Add proteins that are in complexes, bound by small molecules, are in unique molecules, etc.
	# TODO: IMPORTANT


	wisniewskiCounts = validation_data.protein.wisniewski2014Data["avgCounts"]



	plt.figure(figsize = (8.5, 11))

	# maxLine = 1.1 * max(bulkContainer.counts().max(), wisniewskiCounts.max())
	# plt.plot([0, maxLine], [0, maxLine], '--r')
	plt.plot(np.log10(wisniewskiCounts + 1), np.log10(view_validation.counts() + 1), 'o', markeredgecolor = 'k', markerfacecolor = 'none')

	plt.xlabel("log10(Wisniewski 2014 Counts)")
	plt.ylabel("log10(Simulation Average Counts)")
	plt.title("Pearson r: %0.2f" % pearsonr(np.log10(view_validation.counts() + 1), np.log10(wisniewskiCounts + 1))[0])

	# plt.show()

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
