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
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import cPickle
from scipy.stats import pearsonr

import mpld3
from mpld3 import plugins, utils

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

	ids_complexation = sim_data.process.complexation.moleculeNames
	ids_complexation_complexes = sim_data.process.complexation.ids_complexes
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes
	ids_translation = sim_data.process.translation.monomerData["id"].tolist()
	ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))
	bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
	view_complexation = bulkContainer.countsView(ids_complexation)
	view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
	view_equilibrium = bulkContainer.countsView(ids_equilibrium)
	view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
	view_translation = bulkContainer.countsView(ids_translation)
	view_validation = bulkContainer.countsView(validation_data.protein.wisniewski2014Data["monomerId"].tolist())
	view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], np.int)
	proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]
	bulkMolecules.close()

	# Account for monomers
	bulkContainer.countsIs(proteinCountsBulk.mean(axis = 0))

	# Account for unique molecules
	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
	nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
	nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
	uniqueMoleculeCounts.close()
	bulkContainer.countsInc(nActiveRibosome.mean(), sim_data.moleculeGroups.s30_fullComplex + sim_data.moleculeGroups.s50_fullComplex)
	bulkContainer.countsInc(nActiveRnaPoly.mean(), sim_data.moleculeGroups.rnapFull)

	# Account for small-molecule bound complexes
	view_equilibrium.countsInc(
		np.dot(sim_data.process.equilibrium.stoichMatrixMonomers(), view_equilibrium_complexes.counts() * -1)
		)

	# Account for monomers in complexed form
	view_complexation.countsInc(
		np.dot(sim_data.process.complexation.stoichMatrixMonomers(), view_complexation_complexes.counts() * -1)
		)

	wisniewskiCounts = validation_data.protein.wisniewski2014Data["avgCounts"]
	proteinIds = validation_data.protein.wisniewski2014Data["monomerId"].tolist()

	fig, ax = plt.subplots(2, sharey=True, figsize = (8.5, 11))

	# Wisniewski Counts
	points = ax[0].scatter(np.log10(wisniewskiCounts + 1), np.log10(view_validation.counts() + 1), c='w', edgecolor = 'k', alpha=.7)
	ax[0].set_xlabel("log10(Wisniewski 2014 Counts)")
	ax[0].set_title("Pearson r: %0.2f" % pearsonr(np.log10(view_validation.counts() + 1), np.log10(wisniewskiCounts + 1))[0])

	labels = list(proteinIds)
	tooltip = plugins.PointLabelTooltip(points, labels)
	plugins.connect(fig, tooltip)

	# Schmidt Counts
	schmidtLabels = validation_data.protein.schmidt2015Data["monomerId"]
	schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]
	schmidtPoints = ax[1].scatter(
		np.log10(schmidtCounts + 1),
		np.log10(view_validation_schmidt.counts() + 1),
		c='w', edgecolor = 'k', alpha=.7)
	ax[1].set_xlabel("log10(Schmidt 2015 Counts)")
	ax[1].set_title("Pearson r: %0.2f" % pearsonr(np.log10(view_validation_schmidt.counts() + 1), np.log10(schmidtCounts + 1))[0])

	tooltip = plugins.PointLabelTooltip(schmidtPoints, list(schmidtLabels))
	plugins.connect(fig, tooltip)

	plt.ylabel("log10(Simulation Average Counts)")
	# NOTE: This Pearson correlation goes up (at the time of writing) about 0.05 if you only
	# include proteins that you have translational efficiencies for
	plt.xlim(xmin=0)
	plt.ylim(ymin=0)

	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("validationDataFile", help = "File containing loaded and pickled validation data", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
