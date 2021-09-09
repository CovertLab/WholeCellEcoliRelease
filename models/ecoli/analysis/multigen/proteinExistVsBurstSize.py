from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all ids reqiured
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomer_data['id'].size
		n_sims = ap.n_generation

		monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = bool)
		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = bool)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = int)

		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			## READ DATA ##
			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			proteinMonomerCounts = monomerCounts.readColumn("monomerCounts")

			## CALCULATIONS ##
			# Calculate if monomer exists over course of cell cycle
			monomerExist = proteinMonomerCounts.sum(axis=0) > 1

			# Calculate if monomer comes close to doubling
			ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(float) + 1)
			# monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)

			# Load transcription initiation event data
			rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

			# Map transcription initiation events to monomers
			initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.cistron_to_monomer_mapping]

			# Log data
			monomerExistMultigen[gen_idx,:] = monomerExist
			ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer


		existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
		averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen.mean(axis=0)
		averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)

		# import ipdb; ipdb.set_trace()
		# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		# probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		# probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		# for idx, burstSize in enumerate(uniqueBurstSizes):
		# 	mask = initiationEventsPerMonomerMultigen == burstSize
		# 	mask_sum = mask.sum()
		# 	probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
		# 	probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())


		# fig, axesList = plt.subplots(4,1)

		scatterAxis = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		xhistAxis = plt.subplot2grid((4,4), (0,0), colspan=3, sharex = scatterAxis)
		yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

		xhistAxis.xaxis.set_visible(False)
		yhistAxis.yaxis.set_visible(False)

		scatterAxis.scatter(averageInitiationEventsPerMonomer, existFractionPerMonomer, marker = '.', alpha = 0.9, lw = 0.1, s = 5)
		scatterAxis.set_ylabel("Fraction of generations protein exists " + r'$(\frac{\sum n_{exist}}{n_{gen}})$')
		scatterAxis.set_xlabel("Average number of transcription events per monomer " + r"$((\sum_{gen=0}^{n_{gen}} n_{init}) / n_{gen})$")

		scatterAxis.set_ylim([0., 1.1])
		scatterAxis.set_xlim([-10., 1000.])

		yhistAxis.hist(existFractionPerMonomer, bins = 125, orientation='horizontal', range = [0., 1.1], log = True)
		xhistAxis.hist(averageInitiationEventsPerMonomer, bins = 125, log = True, range = [-10., 1000.])
		#xhistAxis.set_xscale("log")


		# axesList[0].semilogy(uniqueBurstSizes, probExistByBurstSize)
		# axesList[1].semilogy(uniqueBurstSizes, probDoubleByBurstSize)

		# # axesList[0].set_ylabel("Probability exists")
		# # axesList[1].set_ylabel("Probability doubles")
		# # axesList[1].set_xlabel("Number of transcription events per generation")

		# axesList[2].semilogy(uniqueBurstSizes, probExistByBurstSize)
		# axesList[2].set_xlim([0., 10.])
		# #axesList[2].set_ylim([0.96, 1.0])
		# axesList[3].semilogy(uniqueBurstSizes, probDoubleByBurstSize)
		# axesList[3].set_xlim([0., 10.])
		# #axesList[3].set_ylim([0.96, 1.0])

		# axesList[0].set_ylabel("Probability\nexists")
		# axesList[1].set_ylabel("Probability\ndoubles")
		# axesList[2].set_ylabel("Probability\nexists")
		# axesList[3].set_ylabel("Probability\ndoubles")
		# axesList[3].set_xlabel("Number of transcription events per generation")

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
