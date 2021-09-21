from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure

from six.moves import cPickle
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

		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = int)

		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			## READ DATA ##
			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			proteinMonomerCounts = monomerCounts.readColumn("monomerCounts")

			## CALCULATIONS ##
			# Calculate if monomer comes close to doubling
			ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(float) + 1)

			# Load transcription initiation event data
			rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			initiation_events_per_cistron = rnapData.readColumn("rna_init_event_per_cistron").sum(axis = 0)

			# Map transcription initiation events to monomers
			initiationEventsPerMonomer = initiation_events_per_cistron[sim_data.relation.cistron_to_monomer_mapping]

			# Log data
			ratioFinalToInitialCount[np.logical_and(ratioFinalToInitialCount == 1., proteinMonomerCounts[0,:] == 0)] = np.nan

			ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer

		uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		ratioByBurstSize = np.zeros(uniqueBurstSizes.size)

		burstSizeToPlot = np.zeros(0)
		ratioToPlot = np.zeros(0)
		for idx, burstSize in enumerate(uniqueBurstSizes):
			mask = initiationEventsPerMonomerMultigen == burstSize
			burstSizeToPlot = np.hstack((burstSizeToPlot, np.ones(mask.sum()) * burstSize))
			ratioToPlot = np.hstack((ratioToPlot, ratioFinalToInitialCountMultigen[mask]))


		real_values_mask = np.logical_not(np.logical_or(np.isnan(ratioToPlot), np.isinf(ratioToPlot)))
		burstSizeToPlot = burstSizeToPlot[real_values_mask]
		ratioToPlot = ratioToPlot[real_values_mask]

		mean = ratioToPlot.mean()
		std = ratioToPlot.std()

		scatterAxis = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		xhistAxis = plt.subplot2grid((4,4), (0,0), colspan=3, sharex = scatterAxis)
		yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

		xhistAxis.xaxis.set_visible(False)
		yhistAxis.yaxis.set_visible(False)

		scatterAxis.semilogx(burstSizeToPlot, ratioToPlot, marker = '.', alpha = 0.75, lw = 0.05) # s = 0.5,
		scatterAxis.set_ylabel("Protein monomer " + r"$\frac{count_f + 1}{count_i + 1}$" + "\n" + r"Clipped $[0, 10]$")
		scatterAxis.set_xlabel("Transcription events per generation\n" + r"Clipped $[0, 1000]$")

		scatterAxis.set_ylim([0., 10.])
		scatterAxis.set_xlim([1., 1000.])

		yhistAxis.hist(ratioToPlot, bins = 125, orientation='horizontal', range = [0., 10.], log = True)
		xhistAxis.hist(burstSizeToPlot, bins = np.logspace(0., 10., 125), log = True)#, range = [0., 1000.])
		xhistAxis.set_xscale("log")

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
