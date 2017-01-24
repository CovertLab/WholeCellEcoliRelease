#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

import cPickle

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
FROM_CACHE = False

CLOSE_TO_DOUBLE = 0.1

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all ids reqiured
	sim_data = cPickle.load(open(simDataFile, "rb"))

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	first_build = True

	# Pre-allocate variables. Rows = Generations, Cols = Monomers
	n_monomers = sim_data.process.translation.monomerData['id'].size
	n_sims = ap.n_generation

	#initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)

	for gen_idx, simDir in enumerate([allDir[0]]):
		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		# Load transcription initiation event data
		rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
		initCounts = rnapData.readColumn("rnaInitEvent")
		initCountsPerMonomer = initCounts[:,sim_data.relation.rnaIndexToMonomerMapping]
		#initCountsPerMonomer.sort(axis = 1)

		#initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

		# Map transcription initiation events to monomers
		#initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rnaIndexToMonomerMapping]

		# Log data
		#initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer

	scatterAxis = plt.subplot2grid((4,3), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
	xhistAxis = plt.subplot2grid((4,3), (0,0), colspan=3, sharex = scatterAxis)
	# yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

	xhistAxis.xaxis.set_visible(False)
	# yhistAxis.yaxis.set_visible(False)

	scatterAxis.imshow(initCountsPerMonomer > 0, cmap='binary')#, interpolation='nearest')
	scatterAxis.set_ylabel("Simulation step")
	scatterAxis.set_xlabel("Monomer")


	scatterAxis.set_xlim([0, initCountsPerMonomer.shape[1]])
	scatterAxis.set_ylim([initCountsPerMonomer.shape[0], 0])

	# yhistAxis.barh(left = range(initCountsPerMonomer.shape[0]), height = initCountsPerMonomer.sum(axis=1))
	xhistAxis.bar(range(initCountsPerMonomer.shape[1]), initCountsPerMonomer.sum(axis=0))
	xhistAxis.set_yscale("log")


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "fixme")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
