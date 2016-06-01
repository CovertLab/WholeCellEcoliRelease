#!/usr/bin/env python
"""
Plot protein monomer counts

@author: John Mason, Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of proteins from the KB

	sim_data = cPickle.load(open(simDataFile, "rb"))

	ids_complexation = sim_data.process.complexation.moleculeNames
	ids_complexation_complexes = [ids_complexation[i] for i in np.where((sim_data.process.complexation.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_equilibrium_complexes = [ids_equilibrium[i] for i in np.where((sim_data.process.equilibrium.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_translation = sim_data.process.translation.monomerData["id"].tolist()
	ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))
	bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
	view_complexation = bulkContainer.countsView(ids_complexation)
	view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
	view_equilibrium = bulkContainer.countsView(ids_equilibrium)
	view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
	view_translation = bulkContainer.countsView(ids_translation)

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

	avgCounts = view_translation.counts()

	relativeCounts = avgCounts / avgCounts.sum()

	expectedCountsArbitrary = normalize(
		sim_data.process.transcription.rnaExpression[sim_data.condition][sim_data.relation.rnaIndexToMonomerMapping] *
		sim_data.process.translation.translationEfficienciesByMonomer /
		(np.log(2) / sim_data.doubling_time.asNumber(units.s) + sim_data.process.translation.monomerData["degRate"].asNumber(1 / units.s))
		)

	expectedCountsRelative = expectedCountsArbitrary / expectedCountsArbitrary.sum()

	plt.figure(figsize = (8.5, 11))

	maxLine = 1.1 * max(np.log10(expectedCountsRelative.max() + 1), np.log10(relativeCounts.max() + 1))
	plt.plot([0, maxLine], [0, maxLine], '--r')
	plt.plot(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1), 'o', markeredgecolor = 'k', markerfacecolor = 'none')

	plt.xlabel("log10(Expected protein distribution (from fitting))")
	plt.ylabel("log10(Actual protein distribution (average over life cycle))")
	plt.title("PCC (of log values): %0.2f" % pearsonr(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1))[0])

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
