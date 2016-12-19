#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)

	max_cells_in_gen = 0
	for genIdx in range(ap.n_generation):
		n_cells = len(ap.get_cells(generation = [genIdx]))
		if n_cells > max_cells_in_gen:
			max_cells_in_gen = n_cells

	fig, axesList = plt.subplots(ap.n_generation, sharex = True)

	growth_rate = np.zeros((max_cells_in_gen, ap.n_generation))

	for genIdx in range(ap.n_generation):
		gen_cells = ap.get_cells(generation = [genIdx])
		for simDir in gen_cells:
			simOutDir = os.path.join(simDir, "simOut")
			growthRate = np.nanmean(TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate"))
			growthRate = (1 / units.s) * growthRate
			growthRate = growthRate.asNumber(1 / units.min)

			growth_rate[np.where(simDir == gen_cells)[0], genIdx] = growthRate

	# Plot initial vs final masses
	if ap.n_generation == 1:
		axesList = [axesList]

	for idx, axes in enumerate(axesList):
		if max_cells_in_gen > 1:
			axes.hist(growth_rate[:,idx].flatten(), np.ceil(np.sqrt(growth_rate[:,idx].size)))
		else:
			axes.plot(growth_rate[:,idx], 1, 'x')
			axes.set_ylim([0, 2])
		axes.axvline(growth_rate[:,idx].mean(), color='k', linestyle='dashed', linewidth=2)
		axes.text(growth_rate[:,idx].mean(), 1, "Mean: %.3f Var: %.3f"%(growth_rate[:,idx].mean(),growth_rate[:,idx].var()))

	axesList[-1].set_xlabel("Growth rate (gDCW/gDCW-min))")
	axesList[ap.n_generation / 2].set_ylabel("Frequency")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

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
