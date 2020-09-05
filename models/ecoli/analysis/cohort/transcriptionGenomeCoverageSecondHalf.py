"""
Plots fraction of mRNAs transcribed (out of all genes to be transcribed) for all seeds.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/29/2016
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt
from six.moves import cPickle, range

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		# Get IDs of mRNAs
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rna_data["id"]
		isMRna = sim_data.process.transcription.rna_data['is_mRNA']
		basalExpression = sim_data.process.transcription.rna_expression["basal"]
		mRnaIds = np.where(isMRna)[0]

		mRnaBasalExpression = np.array([basalExpression[x] for x in mRnaIds])
		mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

		# Sort in order of decreasing basal expression
		descendingOrderIndexing = np.argsort(mRnaBasalExpression)[::-1]
		mRnaBasalExpressionSorted = mRnaBasalExpression[descendingOrderIndexing]
		mRnaNamesSorted = mRnaNames[descendingOrderIndexing]

		# Get all cells in each seed
		ap = AnalysisPaths(variantDir, cohort_plot = True)

		if ap.n_generation == 1:
			print("Only runs for 2 or more cells.")
			return

		second_half_cells = ap.get_cells(generation=range(ap.n_generation//2,ap.n_generation))

		# Get number of mRNAs transcribed
		transcribedFreq = []
		for simDir in second_half_cells:
			simOutDir = os.path.join(simDir, "simOut")

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNamesSorted])
			moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]
			bulkMolecules.close()

			moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
			mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])

			transcribedFreq.append(mRnasTranscribed)

		transcribedFreq = np.array(transcribedFreq)
		transcribedFreqSumOverSeeds = transcribedFreq.sum(axis = 0)

		# Plot
		numMRnas = mRnaNamesSorted.shape[0]
		numCells = second_half_cells.shape[0]
		fig = plt.figure(figsize = (14, 10))

		ax = plt.subplot(1, 1, 1)
		ax.scatter(np.arange(numMRnas), transcribedFreqSumOverSeeds / float(numCells), facecolors = "none", edgecolors = "b")

		heightOffset = 0.01

		ax.set_title("Frequency of producing at least 1 transcript\n(n = %s cells)" % numCells, fontsize = 12)
		ax.set_xlabel("mRNA transcripts\n(in order of decreasing expected basal expression)", fontsize = 10)
		ax.set_xlim([0, numMRnas])
		ax.set_ylim([-.05, 1.05])
		ax.tick_params(which="both", direction="out", top=False)
		ax.spines["top"].set_visible(False)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
