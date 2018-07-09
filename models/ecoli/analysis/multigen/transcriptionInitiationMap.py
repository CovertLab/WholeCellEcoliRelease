from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

FROM_CACHE = False

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
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

		stride = 10
		n_stride = initCountsPerMonomer.shape[0] / stride
		initCountsPerMonomer_downSample = np.zeros((n_stride, n_monomers))
		for idx in range(n_stride):
			initCountsPerMonomer_downSample[idx, :] = initCountsPerMonomer[idx*stride:(idx+1)*stride].sum(axis=0)


		scatterAxis.imshow(initCountsPerMonomer_downSample > 0, cmap='binary', origin = "upper", extent=[0, n_monomers, 0, initCountsPerMonomer_downSample.shape[0]], aspect = initCountsPerMonomer_downSample.shape[0])#, interpolation='nearest')
		scatterAxis.set_ylabel("Simulation step")
		scatterAxis.set_xlabel("Monomer")


		# scatterAxis.set_xlim([0, initCountsPerMonomer_downSample.shape[1]])
		# scatterAxis.set_ylim([initCountsPerMonomer_downSample.shape[0], 0])

		# yhistAxis.barh(left = range(initCountsPerMonomer_downSample.shape[0]), height = initCountsPerMonomer_downSample.sum(axis=1))
		xhistAxis.bar(range(initCountsPerMonomer.shape[1]), initCountsPerMonomer.sum(axis=0))
		xhistAxis.set_yscale("log")

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
