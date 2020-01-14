"""
Plots frequency of observing at least 1 transcript during a cell's life.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/31/2017
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import multigenAnalysisPlot

N_GENES_TO_PLOT = -1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		validation_data = cPickle.load(open(validationDataFile, "rb"))
		essential_RNAs = validation_data.essential_genes.essential_RNAs

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
			simOutDir = os.path.join(simDir, "simOut")

			time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()

			rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			simulatedSynthProb = np.mean(rnaSynthProb.readColumn("rnaSynthProb")[:, mRnaIndexes], axis = 0)
			rnaSynthProb.close()
			simulatedSynthProbs.append(simulatedSynthProb)

			(moleculeCounts,) = read_bulk_molecule_counts(simOutDir, (mRnaIds,))
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
		mRnaIdsOrdered = mRnaIds[indexingOrder]

		alwaysPresentIndexes = np.where(transcribedBoolOrdered == 1.)[0]
		neverPresentIndexes = np.where(transcribedBoolOrdered == 0.)[0]
		sometimesPresentIndexes = np.array([x for x in np.arange(len(transcribedBoolOrdered)) if x not in alwaysPresentIndexes and x not in neverPresentIndexes])
		colors = np.repeat("g", len(transcribedBoolOrdered))
		colors[alwaysPresentIndexes] = "b"
		colors[neverPresentIndexes] = "r"

		# Assemble data
		alwaysTranscriptionEvents_E = []
		alwaysTranscriptionEvents_N = []
		for i in alwaysPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]

			if mRnaIdsOrdered[i] in essential_RNAs:
				alwaysTranscriptionEvents_E.append(v)
			else:
				alwaysTranscriptionEvents_N.append(v)

		neverTranscriptionEvents_E = []
		neverTranscriptionEvents_N = []
		for i in neverPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]

			if mRnaIdsOrdered[i] in essential_RNAs:
				neverTranscriptionEvents_E.append(v)
			else:
				neverTranscriptionEvents_N.append(v)


		sometimesTranscriptionEvents_E = []
		sometimesTranscriptionEvents_N = []
		for i in sometimesPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]

			if mRnaIdsOrdered[i] in essential_RNAs:
				sometimesTranscriptionEvents_E.append(v)
			else:
				sometimesTranscriptionEvents_N.append(v)

		# Plot
		blue = [0, 0, 1]
		green = [0, 0.5, 0]
		red = [1, 0, 0]
		gray = [0, 0, 0]

		fig = plt.figure(figsize = (8, 10))
		scatterAxis = plt.subplot2grid((5, 4), (0, 0), colspan = 3, rowspan = 2)
		histAxis = plt.subplot2grid((5, 4), (0, 3), colspan = 1, rowspan = 2, sharey = scatterAxis)
		alwaysAxis = plt.subplot2grid((5, 4), (2, 0), colspan = 4, rowspan = 1)
		sometimesAxis = plt.subplot2grid((5, 4), (3, 0), colspan = 4, rowspan = 1, sharex = alwaysAxis)
		neverAxis = plt.subplot2grid((5, 4), (4, 0), colspan = 4, rowspan = 1, sharex = alwaysAxis)

		scatterAxis.scatter(np.arange(len(transcribedBoolOrdered)), transcribedBoolOrdered, marker = 'o', facecolors = colors, edgecolors = "none", s = 5)
		scatterAxis.set_title("Frequency of observing at least 1 transcript per generation\n(Genes ordered by simulated synthesis probability)", fontsize = 10)
		scatterAxis.set_xlim([-1, len(transcribedBoolOrdered)])
		scatterAxis.set_ylim([-0.1, 1.1])
		scatterAxis.tick_params(top=False, right=False)
		scatterAxis.tick_params(which = 'both', direction = 'out', labelsize = 8)

		histAxis.hist(transcribedBoolOrdered, bins = len(allDir) + 1, orientation = 'horizontal', color = "k", alpha = 0.5)
		histAxis.set_xscale("log")
		histAxis.spines["right"].set_visible(False)
		histAxis.tick_params(right=False)
		histAxis.tick_params(which = 'both', direction = 'out', labelsize = 8)
		histAxis.text(histAxis.get_xlim()[1] * 1.5, 0, "%s genes\n(%0.1f%%)" % (len(neverTranscriptionEvents_N) + len(neverTranscriptionEvents_E), 100. * (len(neverTranscriptionEvents_N) + len(neverTranscriptionEvents_E)) / float(len(transcribedBoolOrdered))), fontsize = 10, verticalalignment = "top")
		histAxis.text(histAxis.get_xlim()[1] * 1.5, 1, "%s genes\n(%0.1f%%)" % (len(alwaysTranscriptionEvents_N) + len(alwaysTranscriptionEvents_E), 100. * (len(alwaysTranscriptionEvents_N) + len(alwaysTranscriptionEvents_E)) / float(len(transcribedBoolOrdered))), fontsize = 10, verticalalignment = "bottom")
		histAxis.text(histAxis.get_xlim()[1] * 1.5, 0.5, "%s genes\n(%0.1f%%)" % (len(sometimesTranscriptionEvents_N) + len(sometimesTranscriptionEvents_E), 100. * (len(sometimesTranscriptionEvents_N) + len(sometimesTranscriptionEvents_E)) / float(len(transcribedBoolOrdered))), fontsize = 10, verticalalignment = "center")
		histAxis.add_patch(patches.Rectangle((histAxis.get_xlim()[1] * 0.7, 1. / (len(allDir) + 1)), 1e4, 1. - 2. / (len(allDir) + 1), facecolor = green, edgecolor = "none"))

		alwaysAxis.eventplot(alwaysTranscriptionEvents_N + alwaysTranscriptionEvents_E, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [blue] * len(alwaysTranscriptionEvents_N) + [gray] * len(alwaysTranscriptionEvents_E))
		alwaysAxis.set_ylabel("Always present", fontsize = 10)
		alwaysAxis.set_title("Transcription initiation events", fontsize = 10)
		alwaysAxis.set_yticks([])
		alwaysAxis.tick_params(top=False)
		alwaysAxis.tick_params(which = 'both', direction = 'out', labelsize = 8)
		alwaysAxis.set_xlim([0, time[-1] / 3600.])
		alwaysAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(alwaysTranscriptionEvents_E) + len(alwaysTranscriptionEvents_N)])])
		alwaysAxis.text(alwaysAxis.get_xlim()[1] * 1.02, len(alwaysTranscriptionEvents_N) * 0.5, "%s\nnon-essential\ngenes" % len(alwaysTranscriptionEvents_N), fontsize = 10, verticalalignment = "center")
		alwaysAxis.text(alwaysAxis.get_xlim()[1] * 1.02, len(alwaysTranscriptionEvents_N) + len(alwaysTranscriptionEvents_E) * 0.5, "%s essential\ngenes" % len(alwaysTranscriptionEvents_E), fontsize = 10, verticalalignment = "center")

		sometimesAxis.eventplot(sometimesTranscriptionEvents_N + sometimesTranscriptionEvents_E, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [green] * len(sometimesTranscriptionEvents_N) + [gray] * len(sometimesTranscriptionEvents_E))
		sometimesAxis.set_ylabel("Sub-generational", fontsize = 10)
		sometimesAxis.set_yticks([])
		sometimesAxis.tick_params(top=False)
		sometimesAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(sometimesTranscriptionEvents_E) + len(sometimesTranscriptionEvents_N)])])
		sometimesAxis.tick_params(which = 'both', direction = 'out', labelsize = 8)
		sometimesAxis.text(sometimesAxis.get_xlim()[1] * 1.02, len(sometimesTranscriptionEvents_N) * 0.5, "%s\nnon-essential\ngenes" % len(alwaysTranscriptionEvents_N), fontsize = 10, verticalalignment = "center")
		sometimesAxis.text(sometimesAxis.get_xlim()[1] * 1.02, len(sometimesTranscriptionEvents_N) + len(sometimesTranscriptionEvents_E) * 0.5, "%s essential\ngenes" % len(sometimesTranscriptionEvents_E), fontsize = 10, verticalalignment = "center")

		neverAxis.eventplot(neverTranscriptionEvents_N + neverTranscriptionEvents_E, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [red] * len(neverTranscriptionEvents_N) + [gray] * len(neverTranscriptionEvents_E))
		neverAxis.set_ylabel("Never present", fontsize = 10)
		neverAxis.set_xlabel("Time (hour)", fontsize = 10)
		neverAxis.set_yticks([])
		neverAxis.tick_params(top=False)
		neverAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(neverTranscriptionEvents_E) + len(neverTranscriptionEvents_N)])])
		neverAxis.tick_params(which = 'both', direction = 'out', labelsize = 8)
		neverAxis.text(neverAxis.get_xlim()[1] * 1.02, len(neverTranscriptionEvents_N) * 0.5, "%s\nnon-essential\ngenes" % len(neverTranscriptionEvents_N), fontsize = 10, verticalalignment = "center")
		neverAxis.text(neverAxis.get_xlim()[1] * 1.02, len(neverTranscriptionEvents_N) + len(neverTranscriptionEvents_E) * 0.5, "%s essential\ngenes" % len(neverTranscriptionEvents_E), fontsize = 10, verticalalignment = "center")

		plt.subplots_adjust(wspace = 0.4, hspace = 0.4, right = 0.83, bottom = 0.05, left = 0.07, top = 0.95)

		# Only save .png - vectorized formats (.pdf and .svg) are extremely slow
		exportFigure(plt, plotOutDir, plotOutFileName, metadata, extension='.png', dpi=600)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
