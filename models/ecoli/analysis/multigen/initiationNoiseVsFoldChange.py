from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

FROM_CACHE = False

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

CLOSE_TO_DOUBLE = 0.1
FONT_SIZE = 9


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all ids reqiured
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomerData['id'].size
		n_sims = ap.n_generation

		monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = np.bool)
		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = np.float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
		monomerCountInitialMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
		cellMassInitialMultigen = np.zeros(n_sims, dtype = np.float)

		if not FROM_CACHE:

			for gen_idx, simDir in enumerate(allDir):
				simOutDir = os.path.join(simDir, "simOut")

				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				## READ DATA ##
				monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				proteinMonomerCounts = monomerCounts.readColumn("monomerCounts")

				## CALCULATIONS ##
				# Calculate if monomer exists over course of cell cycle
				monomerExist = proteinMonomerCounts.sum(axis=0) > 1

				# Calculate if monomer comes close to doubling
				ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(np.float) + 1)
				# monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)

				# Load transcription initiation event data
				rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
				initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

				# Map transcription initiation events to monomers
				initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rnaIndexToMonomerMapping]

				# Load cell mass
				cellMassInitial = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")[0]

				# Log data
				monomerExistMultigen[gen_idx,:] = monomerExist
				ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
				initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer
				monomerCountInitialMultigen[gen_idx,:] = proteinMonomerCounts[0,:]
				cellMassInitialMultigen[gen_idx] = cellMassInitial

			cPickle.dump(monomerExistMultigen, open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "wb"))
			cPickle.dump(ratioFinalToInitialCountMultigen, open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "wb"))
			cPickle.dump(initiationEventsPerMonomerMultigen, open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "wb"))
			cPickle.dump(monomerCountInitialMultigen, open(os.path.join(plotOutDir,"monomerCountInitialMultigen.pickle"), "wb"))
			cPickle.dump(cellMassInitialMultigen, open(os.path.join(plotOutDir,"cellMassInitialMultigen.pickle"), "wb"))


		monomerExistMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "rb"))
		ratioFinalToInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "rb"))
		initiationEventsPerMonomerMultigen = cPickle.load(open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "rb"))
		monomerCountInitialMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerCountInitialMultigen.pickle"), "rb"))
		cellMassInitialMultigen = cPickle.load(open(os.path.join(plotOutDir,"cellMassInitialMultigen.pickle"), "rb"))
		cellMassInitialMultigen = units.fg * cellMassInitialMultigen

		existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
		averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen.mean(axis=0)
		averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)
		stdInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.std(axis=0)
		cvInitiationEventsPerMonomer = (stdInitiationEventsPerMonomer + 1e-9) / (averageInitiationEventsPerMonomer + 1e-9)

		fig, ax = plt.subplots(1,1)
		ax.semilogy(cvInitiationEventsPerMonomer, ratioFinalToInitialCountMultigen.mean(axis=0), '.', linewidth=0, alpha = 0.25, markeredgewidth=0.)
		ax.set_xlabel("CV of transcription events per monomer " + r"$(\frac{\sigma}{\mu})$", fontsize=FONT_SIZE)
		ax.set_ylabel("Average fold change per\nmonomer over {} generations".format(ap.n_generation), fontsize = FONT_SIZE)
		plt.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			right='off',      # ticks along the bottom edge are off
			left='on',         # ticks along the top edge are off
			)
		whitePadSparklineAxis(ax)

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")
		return

		mws = sim_data.getter.getMass(sim_data.process.translation.monomerData['id'])
		monomerInitialMasses = (mws * monomerCountInitialMultigen / sim_data.constants.nAvogadro)

		# np.tile(cellMassInitialMultigen.asNumber().reshape((1,10)), (n_monomers,1))

		# initialMassFractions = monomerInitialMasses.asNumber(units.fg).transpose() / np.tile(cellMassInitialMultigen.asNumber().reshape((1,10)), (n_monomers,1))
		# averageInitialMassFractions = initialMassFractions.mean(axis = 1)
		# avgMonomerInitialMass = monomerInitialMasses.asNumber(units.fg).mean(axis=0)
		# avgMonomerInitialMassFraction = avgMonomerInitialMass / avgMonomerInitialMass.sum()

		# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		# probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		# probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		# for idx, burstSize in enumerate(uniqueBurstSizes):
		# 	mask = initiationEventsPerMonomerMultigen == burstSize
		# 	mask_sum = mask.sum()
		# 	probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
		# 	probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())


		# fig, axesList = plt.subplots(4,1)

		scatterAxis = plt.subplot2grid((4,5), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		# scatterAxis.axhline(1.0, linewidth=0.5, color='black', linestyle="--", xmin = 0.5, xmax = 1.)
		# scatterAxis.axhline(2.0, linewidth=0.5, color='black', linestyle="--", xmin = 0.5, xmax = 1.)
		xhistAxis = plt.subplot2grid((4,5), (0,0), colspan=3, sharex = scatterAxis)
		yhistAxis = plt.subplot2grid((4,5), (1,3), rowspan=3)#, sharey = scatterAxis)
		yhistAxis.axhline(1.0, linewidth=0.5, color='black', linestyle="--")
		yhistAxis.axhline(2.0, linewidth=0.5, color='black', linestyle="--")
		#yhistAxis_2 = plt.subplot2grid((4,5), (1,4), rowspan=3, sharey = scatterAxis)
		#yhistAxis_2.axhline(1.0, linewidth=0.5, color='black', linestyle="--")
		#yhistAxis_2.axhline(2.0, linewidth=0.5, color='black', linestyle="--")

		xhistAxis.xaxis.set_visible(False)
		yhistAxis.yaxis.set_visible(False)

		smallBurst = averageInitiationEventsPerMonomer <= 1.
		scatterAxis.set_xlim([1e-1, 1e3])
		scatterAxis.set_ylim([0.7, 100])

		# scatterAxis.semilogx(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "red", alpha = 0.9, lw = 0.)#, s = 5)
		# scatterAxis.semilogx(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "blue", alpha = 0.9, lw = 0.)#, s = 5)
		## scatterAxis.semilogx(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "green", alpha = 0.9, lw = 0.)#, s = 5)
		## scatterAxis.semilogx(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "blue", alpha = 0.9, lw = 0.)#, s = 5)

		scatterAxis.loglog(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "green", alpha = 0.9, lw = 0.)#, s = 5)
		scatterAxis.loglog(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "blue", alpha = 0.9, lw = 0.)#, s = 5)


		scatterAxis.set_ylabel("Average fold change per\nmonomer over {} generations".format(ap.n_generation), fontsize = FONT_SIZE)
		scatterAxis.set_xlabel("Average number of transcription\nevents per monomer over {} generations".format(ap.n_generation), fontsize = FONT_SIZE)

		# lims = yhistAxis.get_ylim()
		# step = (lims[1] - lims[0]) / 125
		# bins = np.arange(lims[0], lims[1] + step, step)

		# mass_in_binrange = np.zeros(bins.size-1, dtype=np.float)
		# for i in range(len(bins) - 1):
		# 	in_bin_range = np.logical_and(averageFoldChangePerMonomer > bins[i], averageFoldChangePerMonomer < bins[i+1])
		# 	mass_in_binrange[i] = avgMonomerInitialMassFraction[in_bin_range].sum()

		#yhistAxis_2.barh(bottom = bins[:-1], width = mass_in_binrange, height=(lims[1] - lims[0]) / 125, color = "white")
		#yhistAxis_2.set_xlim([0., 1.])
		#yhistAxis_2.yaxis.set_label_position("right")
		#yhistAxis_2.set_xlabel("Fraction of\nproteome mass", fontsize = FONT_SIZE)
		#scatterAxis.set_xlim([-10., 1000.])

		# yhistAxis.hist(averageFoldChangePerMonomer[~smallBurst], histtype = 'step', bins = 25, orientation='horizontal', log = True)
		# yhistAxis.hist(averageFoldChangePerMonomer[smallBurst], histtype = 'step', bins = 100, orientation='horizontal', log = True, color="green")

		yhistAxis.hist(averageFoldChangePerMonomer[~smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.1), np.log10(100.), 25), range = [0.7, 100], log = True,  orientation='horizontal')
		yhistAxis.hist(averageFoldChangePerMonomer[smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.1), np.log10(100.), 125), range = [0.7, 100], log = True,  orientation='horizontal', color="green")
		yhistAxis.set_ylim([0.7, 100])
		yhistAxis.set_yscale("log")


		xhistAxis.hist(averageInitiationEventsPerMonomer[smallBurst], histtype = 'step', color = "green", bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), log = True, range = [-10., 1000.])
		xhistAxis.hist(averageInitiationEventsPerMonomer[~smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), log = True, range = [-10., 1000.])
		xhistAxis.set_xscale("log")

		for label in yhistAxis.xaxis.get_ticklabels()[::2]:
			label.set_visible(False)

		# for label in yhistAxis_2.xaxis.get_ticklabels()[::2]:
		# 	label.set_visible(False)

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
