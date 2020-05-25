"""
Plots transcription events across multiple generations

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/7/2017
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

N_GENES_TO_PLOT = -1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Get mRNA data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rnaData["id"]
		isMRna = sim_data.process.transcription.rnaData["isMRna"]
		mRnaIndexes = np.where(isMRna)[0]
		mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

		# Get whether or not mRNAs were transcribed
		time = []
		transcribedBool = []
		simulatedSynthProbs = []
		transcriptionEvents = []
		for gen, simDir in enumerate(allDir):
			# print(gen)
			simOutDir = os.path.join(simDir, "simOut")

			time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()

			rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			simulatedSynthProb = np.mean(rnaSynthProb.readColumn("rnaSynthProb")[:, mRnaIndexes], axis = 0)
			rnaSynthProb.close()
			simulatedSynthProbs.append(simulatedSynthProb)

			mRNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'mRNACounts'))
			moleculeCounts = mRNA_counts_reader.readColumn("mRNA_counts")
			moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
			mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])
			transcribedBool.append(mRnasTranscribed)

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			rnaInitEvent = rnapDataReader.readColumn("rnaInitEvent")[:, mRnaIndexes]
			rnapDataReader.close()

			if gen == 0:
				transcriptionEvents = (rnaInitEvent != 0)
			else:
				transcriptionEvents = np.vstack((transcriptionEvents, (rnaInitEvent != 0)))

		time = np.array(time)
		transcribedBool = np.array(transcribedBool)
		simulatedSynthProbs = np.array(simulatedSynthProbs)

		indexingOrder = np.argsort(np.mean(simulatedSynthProbs, axis = 0))
		transcribedBoolOrdered = np.mean(transcribedBool, axis = 0)[indexingOrder]
		transcriptionEventsOrdered = transcriptionEvents[:, indexingOrder]

		alwaysPresentIndexes = np.where(transcribedBoolOrdered == 1.)[0]
		neverPresentIndexes = np.where(transcribedBoolOrdered == 0.)[0]
		sometimesPresentIndexes = np.array([x for x in np.arange(len(transcribedBoolOrdered)) if x not in alwaysPresentIndexes and x not in neverPresentIndexes])
		colors = np.repeat("g", len(transcribedBoolOrdered))
		colors[alwaysPresentIndexes] = "b"
		colors[neverPresentIndexes] = "r"

		# Assemble data
		always=[]
		for i in alwaysPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			always.append(v)

		never=[]
		for i in neverPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			never.append(v)

		sometimes=[]
		for i in sometimesPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			sometimes.append(v)

		# Plot
		blue = [0, 0, 1]
		green = [0, 0.5, 0]

		fig = plt.figure(figsize = (12, 8))
		alwaysAxis = plt.subplot(2, 1, 1)
		sometimesAxis = plt.subplot(2, 1, 2)

		alwaysAxis.eventplot(always, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [blue])
		alwaysAxis.set_ylabel("Frequency == 1", fontsize = 12)
		alwaysAxis.set_xlim([0, time[-1] / 3600.])
		alwaysAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(always)])])
		alwaysAxis.set_xticks([])
		alwaysAxis.set_yticks([])
		alwaysAxis.tick_params(top = False)
		alwaysAxis.tick_params(bottom = False)

		sometimesAxis.eventplot(sometimes, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [green])
		sometimesAxis.set_ylabel("0 < Frequency < 1", fontsize = 12)
		sometimesAxis.set_xlim([0, time[-1] / 3600.])
		sometimesAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(sometimes)])])
		sometimesAxis.set_yticks([])
		sometimesAxis.tick_params(top = False)
		sometimesAxis.tick_params(which = 'both', direction = 'out', labelsize = 12)
		sometimesAxis.set_xticks([0, time[-1] / 3600.])
		sometimesAxis.set_xlabel("Time (hour)", fontsize = 12)

		plt.subplots_adjust(wspace = 0.4, hspace = 0, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)

		# Only save .png - vectorized formats (.pdf and .svg) are extremely slow
		exportFigure(plt, plotOutDir, plotOutFileName, metadata, extension='.png', dpi=600)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
