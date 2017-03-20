#!/usr/bin/env python
"""
Plots Figure 5D.

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
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.utils import units
from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS

BUILD_CACHE = True
def clearLabels(axis):
	axis.set_yticklabels([])
	axis.set_ylabel("")

def bold(lines):
	for line in lines:
		line.set_linewidth(2)

def unbold(lines):
	for line in lines:
		line.set_linewidth(1)

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeComplexId = "MENE-CPLX[c]"
	enzymeMonomerId = "O-SUCCINYLBENZOATE-COA-LIG-MONOMER[c]"
	enzymeRnaId = "EG12437_RNA[c]"
	reactionId = "O-SUCCINYLBENZOATE-COA-LIG-RXN"
	metaboliteIds = ["REDUCED-MENAQUINONE[c]", "CPD-12115[c]"]

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	if 0 not in ap._path_data["seed"]:
		print "Skipping -- figure5D only runs for seed 0"
		return

	allDir = ap.get_cells(seed = [0])

	sim_data = cPickle.load(open(simDataFile, "rb"))
	cellDensity = sim_data.constants.cellDensity
	nAvogadro = sim_data.constants.nAvogadro

	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	mRnaIndexes = np.where(isMRna)[0]
	mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

	simOutDir = os.path.join(allDir[0], "simOut")
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	enzymeComplexIndex = moleculeIds.index(enzymeComplexId)
	enzymeMonomerIndex = moleculeIds.index(enzymeMonomerId)
	enzymeRnaIndex = moleculeIds.index(enzymeRnaId)
	metaboliteIndexes = [moleculeIds.index(x) for x in metaboliteIds]
	bulkMolecules.close()

	if BUILD_CACHE:
		time = []
		enzymeFluxes = []
		enzymeComplexCounts = []
		enzymeMonomerCounts = []
		enzymeRnaCounts = []
		enzymeRnaInitEvent = []
		metaboliteCounts = np.array([])

		cellMass = []
		dryMass = []
		timeStepSec = []

		nTranscriptionInitEventsPerGen = []
		nAvgTetramersPerGen = []

		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()
			timeStepSec += TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec").tolist()
			cellMass += TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass").tolist()
			dryMass += TableReader(os.path.join(simOutDir, "Mass")).readColumn("dryMass").tolist()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")
			enzymeComplexCountsInThisGen = moleculeCounts[:, enzymeComplexIndex].tolist()
			enzymeMonomerCounts += moleculeCounts[:, enzymeMonomerIndex].tolist()
			enzymeRnaCounts += moleculeCounts[:, enzymeRnaIndex].tolist()

			enzymeComplexCounts += enzymeComplexCountsInThisGen
			nAvgTetramersPerGen.append(np.mean(enzymeComplexCountsInThisGen))

			if gen == 0:
				metaboliteCounts = moleculeCounts[:, metaboliteIndexes]
			else:
				metaboliteCounts = np.vstack((metaboliteCounts, moleculeCounts[:, metaboliteIndexes]))
			bulkMolecules.close()

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
			reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
			enzymeFluxes += reactionFluxes[:, np.where(reactionIDs == reactionId)[0][0]].tolist()
			fbaResults.close()

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			rnaInitEventsInThisGen = rnapDataReader.readColumn("rnaInitEvent")[:, np.where(rnaIds == enzymeRnaId)[0][0]].tolist()
			rnapDataReader.close()

			enzymeRnaInitEvent += rnaInitEventsInThisGen
			nTranscriptionInitEventsPerGen.append(np.sum(rnaInitEventsInThisGen))


		time = np.array(time)
		cPickle.dump({
			"time": time,
			"enzymeRnaInitEvent": enzymeRnaInitEvent,
			"enzymeRnaCounts": enzymeRnaCounts,
			"enzymeMonomerCounts": enzymeMonomerCounts,
			"enzymeComplexCounts": enzymeComplexCounts,
			"enzymeFluxes": enzymeFluxes,
			"metaboliteCounts": metaboliteCounts,
			"dryMass": dryMass,
			"cellMass": cellMass,
			"timeStepSec": timeStepSec,
			"nTranscriptionInitEventsPerGen": nTranscriptionInitEventsPerGen,	# storing value to report in paper
			"nAvgTetramersPerGen": nAvgTetramersPerGen,							# storing value to report in paper
			}, open(os.path.join(plotOutDir, "figure5D.pickle"), "wb"))
	else:
		D = cPickle.load(open(os.path.join(plotOutDir, "figure5D.pickle"), "r"))
		time = D["time"]
		enzymeRnaInitEvent = D["enzymeRnaInitEvent"]
		enzymeRnaCounts = D["enzymeRnaCounts"]
		enzymeMonomerCounts = D["enzymeMonomerCounts"]
		enzymeComplexCounts = D["enzymeComplexCounts"]
		enzymeFluxes = D["enzymeFluxes"]
		metaboliteCounts = D["metaboliteCounts"]
		dryMass = D["dryMass"]
		cellMass = D["cellMass"]
		timeStepSec = D["timeStepSec"]

	cellVolume = units.g * np.array(cellMass) / cellDensity
	coefficient = (units.fg * np.array(dryMass)) / (units.fg * np.array(cellMass)) * cellDensity * (timeStepSec * units.s)
	enzymeFluxes = (((COUNTS_UNITS / VOLUME_UNITS) * enzymeFluxes) / coefficient).asNumber(units.mmol / units.g / units.h)

	# Plot
	fig = plt.figure(figsize = (14, 10))
	plt.suptitle("O-succinylbenzoate-CoA ligase downstream behaviors", fontsize = 12)
	rnaInitAxis = plt.subplot(6, 1, 1)
	rnaAxis = plt.subplot(6, 1, 2, sharex = rnaInitAxis)
	monomerAxis = plt.subplot(6, 1, 3, sharex = rnaInitAxis)
	complexAxis = plt.subplot(6, 1, 4, sharex = rnaInitAxis)
	fluxAxis = plt.subplot(6, 1, 5, sharex = rnaInitAxis)
	metAxis = plt.subplot(6, 1, 6)

	rnaInitLine = rnaInitAxis.plot(time / 3600., enzymeRnaInitEvent)
	rnaInitAxis.set_ylabel(r"$menE$" + "\n transcription\nevents", fontsize = 12, rotation = 0)
	rnaInitAxis.yaxis.set_label_coords(-.1, 0.25)
	whitePadSparklineAxis(rnaInitAxis, xAxis = False)

	rnaLine = rnaAxis.plot(time / 3600., enzymeRnaCounts)
	rnaAxis.set_ylabel("menE mRNA\ncounts", fontsize = 12, rotation = 0)
	rnaAxis.yaxis.set_label_coords(-.1, 0.25)
	whitePadSparklineAxis(rnaAxis, xAxis = False)

	monomerLine = monomerAxis.plot(time / 3600., enzymeMonomerCounts)
	monomerAxis.set_ylabel("MenE monomer\ncounts", fontsize = 12, rotation = 0)
	monomerAxis.yaxis.set_label_coords(-.1, 0.25)
	whitePadSparklineAxis(monomerAxis, xAxis = False)

	complexLine = complexAxis.plot(time / 3600., enzymeComplexCounts)
	complexAxis.set_ylabel("MenE tetramer\ncounts", fontsize = 12, rotation = 0)
	complexAxis.yaxis.set_label_coords(-.1, 0.25)
	whitePadSparklineAxis(complexAxis, xAxis = False)

	fluxLine = fluxAxis.plot(time / 3600., enzymeFluxes)
	fluxAxis.set_ylabel("SUCBZL flux\n(mmol/gDCW/hour)", fontsize = 12, rotation = 0)
	fluxAxis.yaxis.set_label_coords(-.1, 0.25)
	whitePadSparklineAxis(fluxAxis, xAxis = False)

	metLine = metAxis.plot(time / 3600., metaboliteCounts[:, 0])
	metAxis.set_ylabel("Menaquinone\ncounts", fontsize = 12, rotation = 0)
	metAxis.yaxis.set_label_coords(-.1, 0.25)
	metAxis.set_xlabel("Time (hour)", fontsize = 12)
	metAxis.set_xlim([0, time[-1] / 3600.])
	whitePadSparklineAxis(metAxis)
	metAxis.set_yticklabels(["%0.1e" % metAxis.get_ylim()[0], "%0.1e" % metAxis.get_ylim()[1]])

	plt.subplots_adjust(hspace = 0.2, right = 0.9, bottom = 0.1, left = 0.15, top = 0.9)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)


	bold(rnaInitLine)
	bold(rnaLine)
	bold(monomerLine)
	bold(complexLine)
	bold(fluxLine)
	bold(metLine)
	exportFigure(plt, plotOutDir, plotOutFileName + "__boldLines", metadata)


	unbold(rnaInitLine)
	unbold(rnaLine)
	unbold(monomerLine)
	unbold(complexLine)
	unbold(fluxLine)
	unbold(metLine)
	plt.suptitle("")
	clearLabels(rnaInitAxis)
	clearLabels(rnaAxis)
	clearLabels(monomerAxis)
	clearLabels(complexAxis)
	clearLabels(fluxAxis)
	clearLabels(metAxis)
	metAxis.set_xticklabels([])
	metAxis.set_xlabel("")
	exportFigure(plt, plotOutDir, plotOutFileName + "__clean", metadata)
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