from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot



class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		# countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')


		# Get all ids reqiured
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomerData['id'].size
		n_sims = ap.n_generation

		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = np.float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)

		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			## READ DATA ##
			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			proteinMonomerCounts = monomerCounts.readColumn("monomerCounts")

			## CALCULATIONS ##
			# Calculate if monomer comes close to doubling
			ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(np.float) + 1)

			# Load transcription initiation event data
			# rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			# initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

			# Map transcription initiation events to monomers
			# initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rnaIndexToMonomerMapping]

			# Log data
			ratioFinalToInitialCount[np.logical_and(ratioFinalToInitialCount == 1., proteinMonomerCounts[0,:] == 0)] = np.nan

			ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
			# initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer

		# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		degradationRates = sim_data.process.transcription.rnaData['degRate'].asNumber(1/units.s)
		degradationRatesByMonomer = degradationRates[sim_data.relation.rnaIndexToMonomerMapping]
		uniqueDegRate = np.unique(degradationRatesByMonomer)

		# burstSizeToPlot = np.zeros(0)
		degRateToPlot = np.zeros(0)
		ratioToPlot = np.zeros(0)
		for idx, degRate in enumerate(uniqueDegRate):
			mask = degradationRatesByMonomer == degRate
			try:
				degRateToPlot = np.hstack((degRateToPlot, np.ones(ratioFinalToInitialCountMultigen[:,mask].size) * degRate))
			except Exception as e:
				# print(e)
				# import ipdb; ipdb.set_trace()
				raise
			ratioToPlot = np.hstack((ratioToPlot, ratioFinalToInitialCountMultigen[:,mask].flatten()))

		real_values_mask = np.logical_not(np.logical_or(np.isnan(ratioToPlot), np.isinf(ratioToPlot)))
		degRateToPlot = degRateToPlot[real_values_mask]
		ratioToPlot = ratioToPlot[real_values_mask]

		mean = ratioToPlot.mean()
		std = ratioToPlot.std()

		scatterAxis = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		xhistAxis = plt.subplot2grid((4,4), (0,0), colspan=3, sharex = scatterAxis)
		yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

		xhistAxis.xaxis.set_visible(False)
		yhistAxis.yaxis.set_visible(False)

		scatterAxis.scatter(degRateToPlot, ratioToPlot, marker = '.', alpha = 0.75, lw = 0.05) # s = 0.5,
		scatterAxis.set_ylabel("Protein monomer " + r"$\frac{count_f + 1}{count_i + 1}$" + "\n" + r"Clipped $[0, 10]$")
		scatterAxis.set_xlabel("Rna degradation rate (1/s)")

		scatterAxis.set_ylim([0., 10.])
		scatterAxis.set_xlim([0., 0.03])

		yhistAxis.hist(ratioToPlot, bins = 125, orientation='horizontal', log = True, range = [0., 10.])
		xhistAxis.hist(degRateToPlot, bins = 125, log = True, range = [0., 0.03])
		#xhistAxis.set_xscale("log")

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
