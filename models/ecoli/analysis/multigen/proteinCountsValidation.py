#!/usr/bin/env python
"""
Compare protein counts to Schmidt 2015 data set

@author: Javier	Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/4/2017
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

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

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
	view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	allDir = ap.get_cells()

	View_Validation_Schmidt = []

	fig = plt.figure(figsize = (4, 4))

	for simDir in allDir:
		print simDir

		simOutDir = os.path.join(simDir, "simOut")

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

		view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())
		View_Validation_Schmidt.append(view_validation_schmidt.counts())

	View_Validation_Schmidt = (np.array(View_Validation_Schmidt)).mean(axis = 0)

	# Schmidt Counts
	schmidtLabels = validation_data.protein.schmidt2015Data["monomerId"]
	schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]

	axis = plt.subplot(1,1,1)

	axis.plot(np.log10(schmidtCounts + 1), np.log10(View_Validation_Schmidt + 1), 'o', color = "black", markersize = 6, alpha = 0.1, zorder = 1, markeredgewidth = 0.0)
	print pearsonr( np.log10(View_Validation_Schmidt + 1), np.log10(schmidtCounts + 1) )[0]

	maxLine = np.ceil(
					max((np.log10(schmidtCounts + 1)).max(),
					(np.log10(View_Validation_Schmidt + 1)).max())
				)
	plt.plot([0, maxLine], [0, maxLine], '-k')

	plt.xlim(xmin=0, xmax=maxLine)
	plt.ylim(ymin=0, ymax=maxLine)

	axis.spines["right"].set_visible(False)
	axis.spines["top"].set_visible(False)
	axis.spines["left"].set_position(("outward", 10))
	axis.spines["bottom"].set_position(("outward", 10))
	axis.tick_params(right = "off")
	axis.tick_params(top = "off")
	axis.tick_params(which = "both", direction = "out")

	axis.set_xlim([-0.07, maxLine])
	axis.set_ylim([-0.07, maxLine])

	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
