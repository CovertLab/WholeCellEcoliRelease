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

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
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

	fig, axesList = plt.subplots(ap.n_generation)

	growth_rate = np.zeros((max_cells_in_gen, ap.n_generation))

	for genIdx in range(ap.n_generation):
		gen_cells = ap.get_cells(generation = [genIdx])
		for simDir in gen_cells:
			simOutDir = os.path.join(simDir, "simOut")
			# growthRate = np.nanmean(TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate"))
			# growthRate = (1 / units.s) * growthRate
			# growthRate = growthRate.asNumber(1 / units.min)

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			doubling_time = (time[-1] - initialTime) * units.s
			averageGrowthRate = np.log(2) / doubling_time.asNumber(units.min)
			growth_rate[np.where(simDir == gen_cells)[0], genIdx] = averageGrowthRate

	
	norm_growth_rate = growth_rate / growth_rate.mean(axis=0)
	norm_max = norm_growth_rate.max()
	norm_min = norm_growth_rate.min()
	nbins = 40
	bin_width = (norm_max - norm_min) / nbins

	# Plot initial vs final masses
	if ap.n_generation == 1:
		axesList = [axesList]

	for idx, axes in enumerate(axesList):
		if max_cells_in_gen > 1:
			axes.hist(norm_growth_rate[:,idx].flatten(), bins=np.arange(norm_growth_rate.min(), norm_growth_rate.max(), bin_width))
		else:
			axes.plot(norm_growth_rate[:,idx], 1, 'x')
			axes.set_ylim([0, 2])

		axes.axvline(norm_growth_rate[:,idx].mean() + norm_growth_rate[:,idx].std(), color='k', linewidth=1, alpha = 0.5)
		axes.axvline(norm_growth_rate[:,idx].mean() - norm_growth_rate[:,idx].std(), color='k', linewidth=1, alpha = 0.5)
		# mean = growth_rate[:,idx].mean()
		# variance = growth_rate[:,idx].var()
		# sd = variance**2.
		# ypos = np.array(axes.get_ylim()).mean()
		#axes.text(growth_rate[:,idx].mean(), ypos, "Mean:{}\nSD:{}\nSD/Mean:{}"%(mean, sd, sd / mean))
		axes.set_xticks([0.5, 0.8, 1.0, 1.2, 1.5])

	axesList[-1].set_xlabel("Normed Growth rate")
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
