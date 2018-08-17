from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure

import cPickle

from wholecell.utils import units
from models.ecoli.analysis import multigenAnalysisPlot

FROM_CACHE = False

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all ids required
		sim_data = cPickle.load(open(simDataFile, "rb"))

		monomerMass = sim_data.getter.getMass(sim_data.process.translation.monomerData['id'])
		nAvogadro = sim_data.constants.nAvogadro

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomerData['id'].size
		n_sims = ap.n_generation

		monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = np.bool)
		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = np.float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
		monomerAvgCountMultigen = np.zeros((n_sims, n_monomers), dtype=np.float)
		cellMassMultigen = np.zeros(n_sims, dtype=np.float)
		monomerInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype=np.int)
		cellMassInitMultigen = np.zeros(n_sims, dtype=np.float)
		cellMassFinalMultigen = np.zeros(n_sims, dtype=np.float)

		if not FROM_CACHE:

			for gen_idx, simDir in enumerate(allDir):
				simOutDir = os.path.join(simDir, "simOut")

				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
				cellMassInit = cellMass[0]
				cellMassFinal = cellMass[-1]
				cellMass = cellMass.mean()

				monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				proteinMonomerCounts = monomerCounts.readColumn("monomerCounts")

				## CALCULATIONS ##
				monomerAverageCount = proteinMonomerCounts.mean(axis=0)

				# Calculate if monomer exists over course of cell cycle
				monomerExist = proteinMonomerCounts.sum(axis=0) > 1

				# Calculate if monomer comes close to doubling
				ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(np.float) + 1)
				# monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)
				monomerInitialCount = proteinMonomerCounts[0,:]

				# Load transcription initiation event data
				rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
				initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

				# Map transcription initiation events to monomers
				initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rnaIndexToMonomerMapping]

				# Log data
				monomerExistMultigen[gen_idx,:] = monomerExist
				ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
				initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer
				monomerAvgCountMultigen[gen_idx, :] = monomerAverageCount
				cellMassMultigen[gen_idx] = cellMass
				monomerInitialCountMultigen[gen_idx,:] = monomerInitialCount

				cellMassInitMultigen[gen_idx] = cellMassInit
				cellMassFinalMultigen[gen_idx] = cellMassFinal

			cPickle.dump(monomerExistMultigen, open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "wb"))
			cPickle.dump(ratioFinalToInitialCountMultigen, open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "wb"))
			cPickle.dump(initiationEventsPerMonomerMultigen, open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "wb"))
			cPickle.dump(monomerAvgCountMultigen, open(os.path.join(plotOutDir,"monomerAvgCountMultigen.pickle"), "wb"))
			cPickle.dump(cellMassMultigen, open(os.path.join(plotOutDir,"cellMassMultigen.pickle"), "wb"))
			cPickle.dump(monomerInitialCountMultigen, open(os.path.join(plotOutDir,"monomerInitialCountMultigen.pickle"), "wb"))
			cPickle.dump(cellMassInitMultigen, open(os.path.join(plotOutDir,"cellMassInitMultigen.pickle"), "wb"))
			cPickle.dump(cellMassFinalMultigen, open(os.path.join(plotOutDir,"cellMassFinalMultigen.pickle"), "wb"))



		monomerExistMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "rb"))
		ratioFinalToInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "rb"))
		initiationEventsPerMonomerMultigen = cPickle.load(open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "rb"))
		monomerAvgCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerAvgCountMultigen.pickle"), "rb"))
		cellMassMultigen = cPickle.load(open(os.path.join(plotOutDir,"cellMassMultigen.pickle"), "rb"))
		monomerInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerInitialCountMultigen.pickle"), "rb"))
		cellMassInitMultigen = cPickle.load(open(os.path.join(plotOutDir,"cellMassInitMultigen.pickle"), "rb"))
		cellMassFinalMultigen = cPickle.load(open(os.path.join(plotOutDir,"cellMassFinalMultigen.pickle"), "rb"))

		monomerFinalCountMultigen = ratioFinalToInitialCountMultigen * (monomerInitialCountMultigen + 1) - 1

		# existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
		# averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen.mean(axis=0)
		# averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)
		# averageCountPerMonomer = monomerAvgCountMultigen.mean(axis=0)

		# translationEff = sim_data.process.translation.translationEfficienciesByMonomer

		initialMonomerTotalMass = ((monomerMass * monomerInitialCountMultigen) / nAvogadro).asNumber(units.fg)
		initialMassFraction = initialMonomerTotalMass / np.tile(cellMassInitMultigen.reshape((ap.n_generation,1)), (1,n_monomers))

		finalMonomerTotalMass = ((monomerMass * monomerFinalCountMultigen) / nAvogadro).asNumber(units.fg)
		finalMassFraction = finalMonomerTotalMass / np.tile(cellMassFinalMultigen.reshape((ap.n_generation,1)), (1,n_monomers))

		massFractionRatio = (finalMassFraction / initialMassFraction).mean(axis=0)
		# Filter out Inf, NaN, and non-positive values before computing log2().
		massFractionRatio_realValues = massFractionRatio[np.isfinite(massFractionRatio) & (massFractionRatio > 0)]
		massFractionRatio_logs = np.log2(massFractionRatio_realValues)

		if massFractionRatio_logs.size > 0:
			fig, axesList = plt.subplots(1)

			axesList.hist(massFractionRatio_logs, bins = 50, log = False)
			axesList.set_xticks(range(15))
			axesList.set_xticklabels(range(15))
			axesList.set_xlabel("Mean mass fraction fold change per monomer")
			axesList.set_ylabel("Count")


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

		# scatterAxis = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		# xhistAxis = plt.subplot2grid((4,4), (0,0), colspan=3, sharex = scatterAxis)
		# yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

		# xhistAxis.xaxis.set_visible(False)
		# yhistAxis.yaxis.set_visible(False)

		# scatterAxis.loglog(averageCountPerMonomer, averageFoldChangePerMonomer, marker = '.', alpha = 0.9, lw = 0.)#, s = 5)
		# scatterAxis.set_ylabel("Average fold change per monomer per generation")
		# scatterAxis.set_xlabel("Average count per generation per monomer")

		# scatterAxis.set_ylim([10e-1, 10e1])
		# scatterAxis.set_xlim([10e-3, 10e5])

		# yhistAxis.hist(averageFoldChangePerMonomer, bins = np.logspace(np.log10(10e-1), np.log10(10e1), 125), orientation='horizontal', log = True, range = [10e-3, 10e6])
		# yhistAxis.set_yscale("log")

		# xhistAxis.hist(averageCountPerMonomer, bins = np.logspace(np.log10(10e-3), np.log10(10e5), 125), log = True, range = [1e-3, 1e5])
		# xhistAxis.set_xscale("log")



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
