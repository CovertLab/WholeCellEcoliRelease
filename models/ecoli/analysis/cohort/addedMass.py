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

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)

	if ap.n_generation == 1:
		print "Need more data to create addedMass"
		return

	max_cells_in_gen = 0
	for genIdx in range(ap.n_generation):
		n_cells = len(ap.get_cells(generation = [genIdx]))
		if n_cells > max_cells_in_gen:
			max_cells_in_gen = n_cells

	# fig, axesList = plt.subplots(ap.n_generation + 1, sharey = True, sharex = True, subplot_kw=dict((("aspect",0.4),("adjustable",'box-forced'))))
	

	fig, axesList = plt.subplots(2,1, sharex = True, sharey = True)
	fig.set_figwidth(5)
	fig.set_figheight(10)

	initial_masses = np.zeros((max_cells_in_gen, ap.n_generation))
	final_masses = np.zeros((max_cells_in_gen, ap.n_generation))

	n_cells = 0
	for genIdx in range(ap.n_generation):
		gen_cells = ap.get_cells(generation = [genIdx])
		if genIdx > 0:
			n_cells += len(gen_cells)
		for simDir in gen_cells:
			try:
				simOutDir = os.path.join(simDir, "simOut")
				mass = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = mass.readColumn("cellMass")

				initial_masses[np.where(simDir == gen_cells)[0], genIdx] = cellMass[0] / 1000.
				final_masses[np.where(simDir == gen_cells)[0], genIdx] = cellMass[-1] / 1000.
				added_masses = final_masses - initial_masses
			except IndexError:
				pass

	# Plot initial vs final masses

	# Plot for all but first generation
	axesList[0].plot(initial_masses[:,1:], added_masses[:,1:], 'x')

	# Plot contours for all but first generation
	H, xedges, yedges = np.histogram2d(initial_masses[:,1:].flatten(), added_masses[:,1:].flatten(), bins=np.round(n_cells/10))
	# counts,ybins,xbins,image = matplotlib.pyplot.hist2d(initial_masses[:,1:].flatten(), added_masses[:,1:].flatten(), bins=20)

	# axesList[1].contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()])
	X, Y = np.meshgrid(xedges, yedges)
	axesList[1].contour(X[:-1,:-1], Y[:-1,:-1], H.transpose())

	axesList[0].set_title("n = {}\nAll generations after first".format(n_cells))
	axesList[0].set_ylabel("Added mass (pg)")

	axesList[1].set_xlabel("Initial mass (pg)")
	axesList[1].set_ylabel("Added mass (pg)")

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
