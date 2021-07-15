from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle, range

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		print("DISABLED")
		return

		# Get all ids reqiured
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomer_data['id'].size
		n_sims = ap.n_generation

		monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = bool)
		monomerDoubleMultigen = np.zeros((n_sims, n_monomers), dtype = bool)
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
			ratioFinalToInitialCount = proteinMonomerCounts[-1:] / proteinMonomerCounts[0,:].astype(float)
			monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)

			# Load transcription initiation event data
			rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

			# Map transcription initiation events to monomers
			initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.RNA_to_monomer_mapping]

			# Log data
			monomerExistMultigen[gen_idx,:] = monomerExist
			monomerDoubleMultigen[gen_idx,:] = monomerDouble
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer

		uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		for idx, burstSize in enumerate(uniqueBurstSizes):
			mask = initiationEventsPerMonomerMultigen == burstSize
			mask_sum = mask.sum()
			probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
			probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())

		# Calculate generational standard deviation (row is generation, col is burst size)
		probExistByBurstSizeGen = np.zeros((ap.n_generation, uniqueBurstSizes.size))
		probDoubleByBurstSizeGen = np.zeros((ap.n_generation, uniqueBurstSizes.size))

		for gen_idx in range(ap.n_generation):
			for idx, burstSize in enumerate(uniqueBurstSizes):
				mask = initiationEventsPerMonomerMultigen[gen_idx,:] == burstSize

				probExistByBurstSize[gen_idx, idx] = monomerExistMultigen[gen_idx,:][mask].sum() / float(mask.sum())

				probDoubleByBurstSize[gen_idx, idx] = monomerDoubleMultigen[gen_idx,:][mask].sum() / float(mask.sum())

		fig, axesList = plt.subplots(4,1)

		axesList[0].semilogy(uniqueBurstSizes, probExistByBurstSize)
		axesList[1].semilogy(uniqueBurstSizes, probDoubleByBurstSize)

		# axesList[0].set_ylabel("Probability exists")
		# axesList[1].set_ylabel("Probability doubles")
		# axesList[1].set_xlabel("Number of transcription events per generation")

		axesList[2].semilogy(uniqueBurstSizes, probExistByBurstSize)
		axesList[2].set_xlim([0., 10.])
		#axesList[2].set_ylim([0.96, 1.0])
		axesList[3].semilogy(uniqueBurstSizes, probDoubleByBurstSize)
		axesList[3].set_xlim([0., 10.])
		#axesList[3].set_ylim([0.96, 1.0])

		axesList[0].set_ylabel("Probability\nexists")
		axesList[1].set_ylabel("Probability\ndoubles")
		axesList[2].set_ylabel("Probability\nexists")
		axesList[3].set_ylabel("Probability\ndoubles")
		axesList[3].set_xlabel("Number of transcription events per generation")


		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")

		# Ignore first 5 generations

		fig, axesList = plt.subplots(4,1)

		uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen[5:, :])
		probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		for idx, burstSize in enumerate(uniqueBurstSizes):
			mask = initiationEventsPerMonomerMultigen[5:, :] == burstSize
			mask_sum = mask.sum()
			probExistByBurstSize[idx] = monomerExistMultigen[5:, :][mask].sum() / float(mask.sum())
			probDoubleByBurstSize[idx] = monomerDoubleMultigen[5:, :][mask].sum() / float(mask.sum())

		axesList[0].plot(uniqueBurstSizes, probExistByBurstSize)
		axesList[1].plot(uniqueBurstSizes, probDoubleByBurstSize)

		# axesList[0].set_ylabel("Probability exists")
		# axesList[1].set_ylabel("Probability doubles")
		# axesList[1].set_xlabel("Number of transcription events per generation")

		axesList[2].plot(uniqueBurstSizes, probExistByBurstSize)
		axesList[2].set_xlim([0., 10.])
		# axesList[2].set_ylim([0.96, 1.0])
		axesList[3].plot(uniqueBurstSizes, probDoubleByBurstSize)
		axesList[3].set_xlim([0., 10.])
		# axesList[3].set_ylim([0.96, 1.0])

		axesList[0].set_ylabel("Probability\nexists")
		axesList[1].set_ylabel("Probability\ndoubles")
		axesList[2].set_ylabel("Probability\nexists")
		axesList[3].set_ylabel("Probability\ndoubles")
		axesList[3].set_xlabel("Number of transcription events per generation")

		exportFigure(plt, plotOutDir, plotOutFileName + "_skip_5_gen",metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
