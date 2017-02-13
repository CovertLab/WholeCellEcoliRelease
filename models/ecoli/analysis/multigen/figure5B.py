#!/usr/bin/env python
"""
Plots Figure 5B.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/12/2017
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

BUILD_CACHE = True

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells(seed = [0])

	# Get mRNA data
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]
	mRnaIndexes = np.where(isMRna)[0]

	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIndexes])
	mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

	if BUILD_CACHE:
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

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			mRnaIndexes_bulk = np.array([moleculeIds.index(x) for x in mRnaIds])
			moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes_bulk]
			bulkMolecules.close()
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
		simulatedSynthProbsOrdered = np.mean(simulatedSynthProbs, axis = 0)[indexingOrder]
		transcriptionEventsOrdered = transcriptionEvents[:, indexingOrder]
		mRnaIdsOrdered = mRnaIds[indexingOrder]

		alwaysPresentIndexes = np.where(transcribedBoolOrdered == 1.)[0]
		neverPresentIndexes = np.where(transcribedBoolOrdered == 0.)[0]
		sometimesPresentIndexes = np.array([x for x in np.arange(len(transcribedBoolOrdered)) if x not in alwaysPresentIndexes and x not in neverPresentIndexes])
		colors = np.repeat("g", len(transcribedBoolOrdered))
		colors[alwaysPresentIndexes] = "b"
		colors[neverPresentIndexes] = "r"
		always = transcribedBoolOrdered[alwaysPresentIndexes]
		sometimes = transcribedBoolOrdered[sometimesPresentIndexes]
		never = transcribedBoolOrdered[neverPresentIndexes]

		alwaysTranscriptionEvents = []
		for i in alwaysPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			alwaysTranscriptionEvents.append(v)
		
		neverTranscriptionEvents = []
		for i in neverPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			neverTranscriptionEvents.append(v)

		sometimesTranscriptionEvents = []
		for i in sometimesPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			sometimesTranscriptionEvents.append(v)

		cPickle.dump({
			"time": time, 
			"transcribedBoolOrdered": transcribedBoolOrdered,
			"colors": colors,
			"id": mRnaIdsOrdered,
			"always": always,
			"sometimes": sometimes,
			"never": never,
			"alwaysTranscriptionEvents": alwaysTranscriptionEvents,
			"sometimesTranscriptionEvents": sometimesTranscriptionEvents,
			"neverTranscriptionEvents": neverTranscriptionEvents,
			}, open(os.path.join(plotOutDir, "figure5B.pickle"), "wb"))

	else:
		D = cPickle.load(open(os.path.join(plotOutDir, "figure5B.pickle"), "r"))
		time = D["time"]
		transcribedBoolOrdered = D["transcribedBoolOrdered"]
		colors = D["colors"]
		mRnaIdsOrdered = D["id"]
		always = D["always"]
		sometimes = D["sometimes"]
		never = D["never"]
		alwaysTranscriptionEvents = D["alwaysTranscriptionEvents"]
		sometimesTranscriptionEvents = D["sometimesTranscriptionEvents"]
		neverTranscriptionEvents = D["neverTranscriptionEvents"]

	# Plot
	fig = plt.figure(figsize = (16, 8))
	plt.suptitle("Frequency of observing at least 1 transcript per generation", fontsize = 12)
	scatterAxis = plt.subplot2grid((2, 4), (0, 0), colspan = 3, rowspan = 2)
	histAxis = plt.subplot2grid((2, 4), (0, 3), colspan = 1, rowspan = 2, sharey = scatterAxis)

	scatterAxis.scatter(np.arange(len(transcribedBoolOrdered)), transcribedBoolOrdered, marker = 'o', facecolors = colors, edgecolors = "none", s = 20)
	scatterAxis.set_xlim([0, len(transcribedBoolOrdered)])
	scatterAxis.set_ylim([0, 1])
	whitePadSparklineAxis(scatterAxis)
	scatterAxis.text(0.5 * scatterAxis.get_xlim()[1], -0.1, "Genes ordered by simulated synthesis probability", horizontalalignment = "center")

	N, bins, patches = histAxis.hist(transcribedBoolOrdered, bins = len(allDir), orientation = 'horizontal')
	for i in xrange(1, len(patches) - 1):
		plt.setp(patches[i], facecolor = "none", edgecolor = "g")
	plt.setp(patches[0], facecolor = "none", edgecolor = "r")
	plt.setp(patches[-1], facecolor = "none", edgecolor = "b")
	histAxis.set_xscale("log")
	whitePadSparklineAxis(histAxis)
	histAxis.xaxis.tick_bottom()

	scatterAxis.set_ylim([-.01, 1.01])

	histAxis.text(histAxis.get_xlim()[1] * 1.7, 0, "%s genes\n(%0.1f%%)" % (len(never), 100. * (len(never) / float(len(transcribedBoolOrdered)))), fontsize = 12, verticalalignment = "center", color = "r")
	histAxis.text(histAxis.get_xlim()[1] * 1.7, 1, "%s genes\n(%0.1f%%)" % (len(always), 100. * (len(always) / float(len(transcribedBoolOrdered)))), fontsize = 12, verticalalignment = "center", color = "b")
	histAxis.text(histAxis.get_xlim()[1] * 1.7, 0.5, "%s genes\n(%0.1f%%)" % (len(sometimes), 100. * (len(sometimes) / float(len(transcribedBoolOrdered)))), fontsize = 12, verticalalignment = "center", color = "g")

	plt.subplots_adjust(wspace = 0.6, hspace = 0.4, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "__top", metadata)
	plt.close("all")


	fig = plt.figure(figsize = (16, 8))
	plt.suptitle("Transcription initiation events", fontsize = 12)
	alwaysAxis = plt.subplot(2, 1, 1)
	sometimesAxis = plt.subplot(2, 1, 2, sharex = alwaysAxis)

	alwaysAxis.eventplot(alwaysTranscriptionEvents, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = "b")
	alwaysAxis.set_ylabel("Frequency == 1", fontsize = 12)
	alwaysAxis.set_xlim([0, time[-1] / 3600.])
	alwaysAxis.set_ylim([-1, len(always)])
	alwaysAxis.set_xticks([])
	alwaysAxis.set_yticks([])
	alwaysAxis.tick_params(top = "off")
	alwaysAxis.tick_params(bottom = "off")
	alwaysAxis.tick_params(axis = "x", labelbottom = 'off')
	
	sometimesAxis.eventplot(sometimesTranscriptionEvents, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = "g")
	sometimesAxis.set_ylabel("0 < Frequency < 1", fontsize = 12)
	sometimesAxis.set_xlim([0, time[-1] / 3600.])
	sometimesAxis.set_ylim([-1, len(sometimes)])
	sometimesAxis.set_xticks([0, time[-1] / 3600.])
	sometimesAxis.set_yticks([])
	sometimesAxis.tick_params(top = "off")
	sometimesAxis.tick_params(which = 'both', direction = 'out', labelsize = 12)
	sometimesAxis.set_xlabel("Time (hour)", fontsize = 12)
	
	plt.subplots_adjust(wspace = 0, hspace = 0, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)
	
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "__bottom", metadata)
	plt.close("all")


if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])