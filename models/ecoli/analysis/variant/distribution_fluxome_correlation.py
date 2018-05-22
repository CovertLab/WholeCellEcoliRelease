#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
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

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

SHUFFLE_VARIANT_TAG = "ShuffleParams"
PLACE_HOLDER = -1

FONT_SIZE=9
trim = 0.05

def getPCC((variant, ap, toyaReactions, toyaFluxesDict, toyaStdevDict)):

	try:

		simDir = ap.get_cells(variant = [variant])[0]

		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
		cellDensity = sim_data.constants.cellDensity

		simOutDir = os.path.join(simDir, "simOut")

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass")
		dryMass = massListener.readColumn("dryMass")
		massListener.close()

		coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
		reactionFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fbaResults.readColumn("reactionFluxes").T / coefficient).T
		fbaResults.close()

		modelFluxes = {}
		for toyaReaction in toyaReactions:
			fluxTimeCourse = []

			for rxn in reactionIDs:
				if re.findall(toyaReaction, rxn):
					reverse = 1
					if re.findall("(reverse)", rxn):
						reverse = -1

					if len(fluxTimeCourse):
						fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
					else:
						fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

			if len(fluxTimeCourse):
				if toyaReaction not in modelFluxes:
					modelFluxes[toyaReaction] = []
				modelFluxes[toyaReaction].append(np.mean(fluxTimeCourse).asNumber(units.mmol / units.g / units.h))

		toyaVsReactionAve = []
		for rxn, toyaFlux in toyaFluxesDict.iteritems():
			if rxn in ["ISOCITDEH-RXN", "SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31."]:
				continue
			if rxn in modelFluxes:
				toyaVsReactionAve.append((np.mean(modelFluxes[rxn]), toyaFlux.asNumber(units.mmol / units.g / units.h), np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(units.mmol / units.g / units.h)))

		toyaVsReactionAve = np.array(toyaVsReactionAve)
		pcc, pval = pearsonr(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])

		return pcc, pval


	except:
		return np.nan, np.nan




def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):

	if metadata is not None and SHUFFLE_VARIANT_TAG not in metadata["variant"]:
		print "This plot only runs for variants where parameters are shuffled."
		return

	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
	toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
	toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
	toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

	ap = AnalysisPaths(inputDir, variant_plot = True)



	pool = Pool(processes = 16)
	args = zip(range(ap.n_variant), [ap] * ap.n_variant, [toyaReactions] * ap.n_variant, [toyaFluxesDict] * ap.n_variant, [toyaStdevDict] * ap.n_variant)
	import time
	start = time.time()
	result = pool.map(getPCC, args)
	end = time.time()
	print end - start
	cPickle.dump(result, open("pcc_results_fluxome.cPickle", "w"), cPickle.HIGHEST_PROTOCOL)
	pool.close()
	pool.join()
	result = cPickle.load(open("pcc_results_fluxome.cPickle", "r"))
	controlPcc, controlPvalue = result[0]
	pccs, pvals = zip(*result[1:])
	pccs = np.array(pccs)
	pvals = np.array(pvals)

	fig = plt.figure()
	fig.set_figwidth(5)
	fig.set_figheight(5)
	ax = plt.subplot(1, 1, 1)

	pccs = np.array([x for x in pccs if not np.isnan(x)])
	ax.hist(pccs, np.sqrt(pccs.size))
	ax.axvline(controlPcc, color = "k", linestyle = "dashed", linewidth = 2)

	ax.set_xlabel("Fluxome correlation (Pearson r)")
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
