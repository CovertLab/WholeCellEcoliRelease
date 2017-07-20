#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as patches


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from scipy.stats import pearsonr
from multiprocessing import Pool

SHUFFLE_VARIANT_TAG = "ShuffleParams"
PLACE_HOLDER = -1

FONT_SIZE=9
trim = 0.05

def getPCC((variant, ap, monomerIds, schmidtCounts)):

	simDir = ap.get_cells(variant = [variant])[0]

	sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))

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
	view_validation_schmidt = bulkContainer.countsView(monomerIds)

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

	pcc, pval = pearsonr(np.log10(view_validation_schmidt.counts() + 1), np.log10(schmidtCounts + 1))

	return pcc, pval




def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):

	if metadata is not None and SHUFFLE_VARIANT_TAG not in metadata["variant"]:
		print "This plot only runs for variants where parameters are shuffled."
		return

	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]

	ap = AnalysisPaths(inputDir, variant_plot = True)


	pool = Pool(processes = 16)
	args = zip(range(ap.n_variant), [ap] * ap.n_variant, [validation_data.protein.schmidt2015Data["monomerId"].tolist()] * ap.n_variant, [schmidtCounts] * ap.n_variant)
	# import time
	# start = time.time()
	result = pool.map(getPCC, args)
	# end = time.time()
	# print end - start
	# cPickle.dump(result, open("pcc_results.cPickle", "w"), cPickle.HIGHEST_PROTOCOL)
	pool.close()
	pool.join()
	# result = cPickle.load(open("pcc_results.cPickle", "r"))
	controlPcc, controlPvalue = result[0]
	pccs, pvals = zip(*result[1:])
	pccs = np.array(pccs)
	pvals = np.array(pvals)

	fig = plt.figure()
	fig.set_figwidth(5)
	fig.set_figheight(5)
	ax = plt.subplot(1, 1, 1)

	ax.hist(pccs, np.sqrt(pccs.size))
	ax.axvline(controlPcc, color = "k", linestyle = "dashed", linewidth = 2)

	ax.set_xlabel("Proteome correlation (Pearson r)")
	ax.set_title("Mean: %0.3g     Std: %0.3g     Control: %0.3g" % (pccs.mean(), pccs.std(), controlPcc))

	axes_list = [ax]

	for a in axes_list:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in a.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 

	whitePadSparklineAxis(ax)

	plt.subplots_adjust(bottom = 0.2, wspace=0.3)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("validationDataFile", help = "Validation file name", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["validationDataFile"])