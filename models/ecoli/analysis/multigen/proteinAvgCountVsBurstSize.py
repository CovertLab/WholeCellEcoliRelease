from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all ids required
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get all cells
		allDir = self.ap.get_cells()

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomer_data['id'].size
		n_sims = self.ap.n_generation

		monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = bool)
		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = int)
		monomerAvgCountMultigen = np.zeros((n_sims, n_monomers), dtype=float)

		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			proteinMonomerCounts = monomerCounts.readColumn("monomerCounts")

			## CALCULATIONS ##
			monomerAverageCount = proteinMonomerCounts.mean(axis=0)

			# Calculate if monomer exists over course of cell cycle
			monomerExist = proteinMonomerCounts.sum(axis=0) > 1

			# Calculate if monomer comes close to doubling
			ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(float) + 1)
			# monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)

			# Load transcription initiation event data
			rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			initiation_events_per_cistron = rnapData.readColumn("rna_init_event_per_cistron").sum(axis = 0)

			# Map transcription initiation events to monomers
			initiationEventsPerMonomer = initiation_events_per_cistron[sim_data.relation.cistron_to_monomer_mapping]

			# Log data
			monomerExistMultigen[gen_idx,:] = monomerExist
			ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer
			monomerAvgCountMultigen[gen_idx, :] = monomerAverageCount

		existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
		averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen.mean(axis=0)
		averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)
		averageCountPerMonomer = monomerAvgCountMultigen.mean(axis=0)

		translationEff = sim_data.process.translation.translation_efficiencies_by_monomer

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

		scatterAxis.loglog(averageInitiationEventsPerMonomer, averageCountPerMonomer, marker = '.', alpha = 0.9, lw = 0.)#, s = 5)
		scatterAxis.set_ylabel("Average count per generation per monomer")
		scatterAxis.set_xlabel("Average number of transcription events per generation per monomer")

		scatterAxis.set_ylim([10e-3, 10e6])
		scatterAxis.set_xlim([10e-3, 1000.])

		yhistAxis.hist(averageCountPerMonomer, bins = np.logspace(np.log10(10e-3), np.log10(10e6), 125), orientation='horizontal', log = True, range = [10e-3, 10e6])
		yhistAxis.set_yscale("log")

		xhistAxis.hist(averageInitiationEventsPerMonomer, bins = np.logspace(np.log10(10e-3), np.log10(10e6), 125), log = True, range = [-10., 1000.])
		xhistAxis.set_xscale("log")



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
