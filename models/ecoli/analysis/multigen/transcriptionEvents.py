#!/usr/bin/env python
"""
Plots transcription events across multiple generations
@author: Sam Bray
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/7/2017
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

USE_CACHE = False
N_GENES_TO_PLOT = -1

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	essentialRnas = validation_data.essentialGenes.essentialRnas

	# Get mRNA data
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]
	mRnaIndexes = np.where(isMRna)[0]

	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIndexes])
	mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

	if not USE_CACHE:
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

		# Assemble data
		alwaysTranscriptionEvents_E = []
		alwaysTranscriptionEvents_N = []
		alwaysId_E = []
		alwaysId_N = []
		always=[]
		for i in alwaysPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			always.append(v)
		
		neverTranscriptionEvents_E = []
		neverTranscriptionEvents_N = []
		never=[]
		for i in neverPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			never.append(v)


		sometimesTranscriptionEvents_E = []
		sometimesTranscriptionEvents_N = []
		sometimes=[]
		for i in sometimesPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			sometimes.append(v)

		cPickle.dump({
			"time": time, 
			"always_E": alwaysTranscriptionEvents_E,
			"always_N": alwaysTranscriptionEvents_N, 
			"never_E": neverTranscriptionEvents_E, 
			"never_N": neverTranscriptionEvents_N,
			"sometimes_E": sometimesTranscriptionEvents_E,
			"sometimes_N": sometimesTranscriptionEvents_N,
			"transcriptionFrequency": transcribedBoolOrdered,
			"colors": colors,
			"id": mRnaIdsOrdered,
			"always": always,
			"never": never,
			"sometimes": sometimes,
			}, open(os.path.join(plotOutDir, "transcriptionEvents.pickle"), "wb"))

	if USE_CACHE:
		D = cPickle.load(open(os.path.join(plotOutDir, "transcriptionEvents.pickle"), "r"))
		time = D["time"]
		alwaysTranscriptionEvents_E = D["always_E"]
		alwaysTranscriptionEvents_N = D["always_N"]
		neverTranscriptionEvents_E = D["never_E"]
		neverTranscriptionEvents_N = D["never_N"]
		sometimesTranscriptionEvents_E = D["sometimes_E"]
		sometimesTranscriptionEvents_N = D["sometimes_N"]
		transcribedBoolOrdered = D["transcriptionFrequency"]
		colors = D["colors"]
		mRnaIdsOrdered = D["id"]
		always = D["always"]
		never = D["never"]
		sometimes = D["sometimes"]

	# Plot
	blue = [0, 0, 1]
	green = [0, 0.5, 0]
	red = [1, 0, 0]
	gray = [0, 0, 0]

	fig = plt.figure(figsize = (8, 10))
	alwaysAxis = plt.subplot2grid((5, 4), (2, 0), colspan = 4, rowspan = 1)
	sometimesAxis = plt.subplot2grid((5, 4), (3, 0), colspan = 4, rowspan = 1, sharex = alwaysAxis)
	neverAxis = plt.subplot2grid((5, 4), (4, 0), colspan = 4, rowspan = 1, sharex = alwaysAxis)

	
	alwaysAxis.eventplot(always, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [blue])
	alwaysAxis.set_title("Transcription initiation events", fontsize = 10)
	alwaysAxis.set_yticks([])
	alwaysAxis.tick_params(top = "off")
	alwaysAxis.tick_params(which = 'both', direction = 'in', labelsize = 0)
	alwaysAxis.set_xlim([0, time[-1] / 3600.])
	alwaysAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(always)])])
	
	sometimesAxis.eventplot(sometimes, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [green])# * len(sometimesTranscriptionEvents_N) + [gray] * len(sometimesTranscriptionEvents_E))
	sometimesAxis.set_ylabel("Sub-generational", fontsize = 10)
	sometimesAxis.set_yticks([])
	sometimesAxis.tick_params(top = "off")
	sometimesAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(sometimes)])])
	sometimesAxis.tick_params(which = 'both', direction = 'in', labelsize = 0)
	
	neverAxis.eventplot(never, orientation = "horizontal", linewidths = 2., linelengths = 1., colors = [red])
	neverAxis.set_ylabel("Never present", fontsize = 10)
	neverAxis.set_xlabel("Time (hour)", fontsize = 10)
	neverAxis.set_yticks([])
	neverAxis.tick_params(top = "off")
	neverAxis.set_ylim([-1, np.max([N_GENES_TO_PLOT, len(never)])])
	neverAxis.tick_params(which = 'both', direction = 'out', labelsize = 8)
	
	plt.subplots_adjust(wspace = 0.4, hspace = 0, right = 0.83, bottom = 0.05, left = 0.07, top = 0.95)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
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
