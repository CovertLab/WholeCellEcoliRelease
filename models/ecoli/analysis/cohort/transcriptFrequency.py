"""
Plots transcript frequency (ie. frequency of observing at least 
one copy of transcript) at the 4th generation across 32 seeds.
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

N_SEEDS = 32


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		# Get all cells in each seed
		allDir = self.ap.get_cells(generation = [3])

		if len(allDir) <= 1:
			print("Skipping -- transcriptFrequency only runs for multiple seeds")
			return

		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get mRNA data
		rnaIds = sim_data.process.transcription.rna_data["id"]
		isMRna = sim_data.process.transcription.rna_data['is_mRNA']
		mRnaIndexes = np.where(isMRna)[0]
		synthProb = sim_data.process.transcription.rna_synth_prob["basal"]

		mRnaSynthProb = np.array([synthProb[x] for x in mRnaIndexes])
		mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

		# Get mRNA indices in bulk molecules
		simOutDir = os.path.join(allDir[0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mRnaIndexes_bulk = [moleculeIds.index(x) for x in mRnaIds]

		# Get frequency data
		hadTranscribed = None
		simulatedSynthProb = None
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes_bulk]
			moleculeCounts_sumOverTime = moleculeCounts.sum(axis = 0)
			mRnasTranscribed = np.array([x != 0 for x in moleculeCounts_sumOverTime])

			rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			simulatedSynthProb_ = rnaSynthProb.readColumn("rnaSynthProb")

			if hadTranscribed is None:
				hadTranscribed = mRnasTranscribed
				simulatedSynthProb = np.mean(simulatedSynthProb_, axis = 0)[isMRna]
				continue
			hadTranscribed = np.vstack((hadTranscribed, mRnasTranscribed))
			simulatedSynthProb = np.vstack((simulatedSynthProb, np.mean(simulatedSynthProb_, axis = 0)[isMRna]))

			if hadTranscribed.shape[0] == N_SEEDS:
				break

		hadTranscribedFrequency = np.mean(hadTranscribed, axis = 0)
		simulatedSynthProbAvg = np.mean(simulatedSynthProb, axis = 0)
		indexOrder = np.argsort(simulatedSynthProbAvg)

		# Plot
		fig = plt.figure(figsize = (16, 8))
		scatterAxis = plt.subplot2grid((2, 4), (0, 0), colspan = 3, rowspan = 2)
		histAxis = plt.subplot2grid((2, 4), (0, 3), colspan = 1, rowspan = 2, sharey = scatterAxis)

		colors = np.repeat("g", len(hadTranscribedFrequency))
		colors[np.where(hadTranscribedFrequency[indexOrder] == 1)] = "b"
		colors[np.where(hadTranscribedFrequency[indexOrder] == 0)] = "r"

		scatterAxis.scatter(np.arange(len(hadTranscribedFrequency)), hadTranscribedFrequency[indexOrder], marker = "o", facecolors = colors, edgecolors = "none", s = 20)
		scatterAxis.set_xlim([0, len(hadTranscribedFrequency)])
		scatterAxis.set_ylim([0, 1])
		whitePadSparklineAxis(scatterAxis)

		N, bins, patches = histAxis.hist(hadTranscribedFrequency, bins = N_SEEDS + 1, orientation = 'horizontal')
		for i in range(1, len(patches) - 1):
			plt.setp(patches[i], facecolor = "none", edgecolor = "g")
		plt.setp(patches[0], facecolor = "none", edgecolor = "r")
		plt.setp(patches[-1], facecolor = "none", edgecolor = "b")
		histAxis.set_xscale("log")
		whitePadSparklineAxis(histAxis)
		histAxis.xaxis.tick_bottom()
		histXmin, histXmax = histAxis.get_xlim()

		scatterAxis.set_ylim([-.01, 1.01])
		plt.subplots_adjust(wspace = 0.6, hspace = 0.4, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close()


if __name__ == "__main__":
	Plot().cli()
