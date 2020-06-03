"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/12/2017
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.utils import units
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from models.ecoli.analysis import multigenAnalysisPlot

FONTSIZE = 6
LABELSIZE = 6
PLOT_DOWNSTREAM = True

def clearLabels(axis):
	axis.set_yticklabels([])
	axis.set_ylabel("")

def bold(lines):
	for line in lines:
		line.set_linewidth(2)

def unbold(lines):
	for line in lines:
		line.set_linewidth(1)

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		enzymeComplexId = "MENE-CPLX[c]"
		enzymeMonomerId = "O-SUCCINYLBENZOATE-COA-LIG-MONOMER[c]"
		enzymeRnaId = "EG12437_RNA[c]"
		reactionId = "O-SUCCINYLBENZOATE-COA-LIG-RXN"
		metaboliteIds = ["REDUCED-MENAQUINONE[c]", "CPD-12115[c]"]

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		if 0 not in ap._path_data["seed"]:
			print("Skipping -- figure5D only runs for seed 0")
			return

		allDir = ap.get_cells(seed = [0])

		sim_data = cPickle.load(open(simDataFile, "rb"))
		cellDensity = sim_data.constants.cellDensity

		rnaIds = sim_data.process.transcription.rnaData["id"]

		simOutDir = os.path.join(allDir[0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		enzymeComplexIndex = moleculeIds.index(enzymeComplexId)
		enzymeMonomerIndex = moleculeIds.index(enzymeMonomerId)
		metaboliteIndexes = [moleculeIds.index(x) for x in metaboliteIds]

		mRNA_counts_reader = TableReader(
			os.path.join(simOutDir, 'mRNACounts'))
		all_mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')
		enzymeRnaIndex = all_mRNA_ids.index(enzymeRnaId)

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
		generationTicks = [0.]

		nTranscriptionInitEventsPerGen = []
		nAvgTetramersPerGen = []

		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			time += main_reader.readColumn("time").tolist()
			generationTicks.append(time[-1])
			timeStepSec += main_reader.readColumn("timeStepSec").tolist()

			mass = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass += mass.readColumn("cellMass").tolist()
			dryMass += mass.readColumn("dryMass").tolist()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")
			enzymeComplexCountsInThisGen = moleculeCounts[:, enzymeComplexIndex].tolist()
			enzymeMonomerCounts += moleculeCounts[:, enzymeMonomerIndex].tolist()
			enzymeComplexCounts += enzymeComplexCountsInThisGen
			nAvgTetramersPerGen.append(np.mean(enzymeComplexCountsInThisGen))

			mRNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'mRNACounts'))
			mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
			enzymeRnaCounts += mRNA_counts[:, enzymeRnaIndex].tolist()

			if gen == 0:
				metaboliteCounts = moleculeCounts[:, metaboliteIndexes]
			else:
				metaboliteCounts = np.vstack((metaboliteCounts, moleculeCounts[:, metaboliteIndexes]))

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
			reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
			enzymeFluxes += reactionFluxes[:, np.where(reactionIDs == reactionId)[0][0]].tolist()

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			rnaInitEventsInThisGen = rnapDataReader.readColumn("rnaInitEvent")[:, np.where(rnaIds == enzymeRnaId)[0][0]].tolist()

			enzymeRnaInitEvent += rnaInitEventsInThisGen
			nTranscriptionInitEventsPerGen.append(np.sum(rnaInitEventsInThisGen))

		time = np.array(time)

		coefficient = (units.fg * np.array(dryMass)) / (units.fg * np.array(cellMass)) * cellDensity * (timeStepSec * units.s)
		enzymeFluxes = (((COUNTS_UNITS // VOLUME_UNITS) * enzymeFluxes) / coefficient).asNumber(units.mmol / units.g / units.h)

		averages = []
		indices = [np.where(time == x)[0][0] for x in generationTicks]
		for x in np.arange(len(indices) - 1):
			avg = np.average(enzymeComplexCounts[indices[x]:indices[x+1]])
			averages.append(avg)

		# Plot
		fig = plt.figure(figsize = (7, 7))
		plt.suptitle("O-succinylbenzoate-CoA ligase downstream behaviors", fontsize = FONTSIZE
			)
		rnaInitAxis = plt.subplot(6, 1, 1)
		rnaAxis = plt.subplot(6, 1, 2, sharex = rnaInitAxis)
		monomerAxis = plt.subplot(6, 1, 3, sharex = rnaInitAxis)
		complexAxis = plt.subplot(6, 1, 4, sharex = rnaInitAxis)
		fluxAxis = plt.subplot(6, 1, 5, sharex = rnaInitAxis)
		metAxis = plt.subplot(6, 1, 6)

		rnaInitLine = rnaInitAxis.plot(time / 3600., enzymeRnaInitEvent)
		rnaInitAxis.set_ylabel(r"$menE$" + "\n transcription\nevents", fontsize = FONTSIZE, rotation = 0)
		rnaInitAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(rnaInitAxis, xAxis = False)

		rnaLine = rnaAxis.plot(time / 3600., enzymeRnaCounts)
		rnaAxis.set_ylabel("menE mRNA\ncounts", fontsize = FONTSIZE, rotation = 0)
		rnaAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(rnaAxis, xAxis = False)

		monomerLine = monomerAxis.plot(time / 3600., enzymeMonomerCounts)
		monomerAxis.set_ylabel("MenE monomer\ncounts", fontsize = FONTSIZE, rotation = 0)
		monomerAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(monomerAxis, xAxis = False)
		monomerAxis.set_yticks([0, 4, monomerAxis.get_ylim()[1]])
		monomerAxis.set_yticklabels(["0", "4", "%s" % monomerAxis.get_ylim()[1]])

		complexLine = complexAxis.plot(time / 3600., enzymeComplexCounts)
		complexAxis.set_ylabel("MenE tetramer\ncounts", fontsize = FONTSIZE, rotation = 0)
		complexAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(complexAxis, xAxis = False)

		fluxLine = fluxAxis.plot(time / 3600., enzymeFluxes)
		fluxAxis.set_ylabel("SUCBZL flux\n(mmol/gDCW/hour)", fontsize = FONTSIZE, rotation = 0)
		fluxAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(fluxAxis, xAxis = False)

		metLine = metAxis.plot(time / 3600., np.sum(metaboliteCounts, axis = 1))
		metAxis.set_ylabel("End product\ncounts", fontsize = FONTSIZE, rotation = 0)
		metAxis.yaxis.set_label_coords(-.1, 0.25)
		metAxis.set_xlabel("Time (hour)\ntickmarks at each new generation", fontsize = FONTSIZE)
		metAxis.set_ylim([metAxis.get_ylim()[0] * 0.2, metAxis.get_ylim()[1]])
		whitePadSparklineAxis(metAxis)
		metAxis.set_yticklabels(["%0.1e" % metAxis.get_ylim()[0], "%0.1e" % metAxis.get_ylim()[1]])
		metAxis.set_xticks(np.array(generationTicks) / 3600.)
		xticklabels = np.repeat("     ", len(generationTicks))
		xticklabels[0] = "0"
		xticklabels[-1] = "%0.2f" % (time[-1] / 3600.)
		metAxis.set_xticklabels(xticklabels)

		noComplexIndexes = np.where(np.array(enzymeComplexCounts) == 0)[0]
		patchStart = []
		patchEnd = []
		if len(noComplexIndexes):
			prev = noComplexIndexes[0]
			patchStart.append(prev)
			for i in noComplexIndexes:
				if np.abs(i - prev) > 1:
					patchStart.append(i)
					patchEnd.append(prev)
				prev = i
			patchEnd.append(prev)

		axesList = [rnaInitAxis, rnaAxis, monomerAxis, complexAxis, fluxAxis, metAxis]
		for axis in axesList:
			axis.tick_params(labelsize = LABELSIZE)
			for i in xrange(len(patchStart)):
				width = time[patchEnd[i]] / 3600. - time[patchStart[i]] / 3600.
				if width <= 0.1:
					continue

				height = axis.get_ylim()[1] - axis.get_ylim()[0]
				axis.add_patch(patches.Rectangle((time[patchStart[i]] / 3600., axis.get_ylim()[0]), width, height, alpha = 0.25, color = "gray", linewidth = 0.))

		plt.subplots_adjust(hspace = 0.5, right = 0.9, bottom = 0.1, left = 0.15, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		lineList = [rnaInitLine, rnaLine, monomerLine, complexLine, fluxLine, metLine]
		axesList = [rnaInitAxis, rnaAxis, monomerAxis, complexAxis, fluxAxis, metAxis]
		for l in lineList:
			bold(l)
		exportFigure(plt, plotOutDir, plotOutFileName + "__boldLines", metadata)

		for l in lineList:
			unbold(l)
		for a in axesList:
			clearLabels(a)
		plt.suptitle("")
		metAxis.set_xticklabels([])
		metAxis.set_xlabel("")
		exportFigure(plt, plotOutDir, plotOutFileName + "__clean", "")
		plt.close("all")

		if PLOT_DOWNSTREAM:
			fig, axesList = plt.subplots(12, figsize = (14, 14))
			plt.subplots_adjust(hspace = 0.5, right = 0.95, bottom = 0.05, left = 0.15, top = 0.95)
			enzymeIds = ["MENE-CPLX[c]", "CPLX0-7882[c]", "CPLX0-8128[c]", "DMK-MONOMER[i]", "2-OCTAPRENYL-METHOXY-BENZOQ-METH-MONOMER[c]"]
			reactionIds = ["O-SUCCINYLBENZOATE-COA-LIG-RXN", "NAPHTHOATE-SYN-RXN", "RXN-9311", "DMK-RXN", "ADOMET-DMK-METHYLTRANSFER-RXN"]
			reactantIds = ["CPD-12115[c]"]

			for gen, simDir in enumerate(allDir):
				simOutDir = os.path.join(simDir, "simOut")

				main_reader = TableReader(os.path.join(simOutDir, "Main"))
				time_ = main_reader.readColumn("time")
				timeStepSec = main_reader.readColumn("timeStepSec")

				mass = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = mass.readColumn("cellMass")
				dryMass = mass.readColumn("dryMass")

				(enzymeCounts, metCounts, reactantCounts) = read_bulk_molecule_counts(
					simOutDir, (enzymeIds, metaboliteIds[:1], reactantIds))

				fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
				reactionIDs_ = np.array(fbaResults.readAttribute("reactionIDs")).tolist()  # type: list
				reactionIndexes = [reactionIDs_.index(x) for x in reactionIds]
				reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
				enzymeFluxes = reactionFluxes[:, reactionIndexes]
				fbaResults.close()

				coefficient = (units.fg * np.array(dryMass)) / (units.fg * np.array(cellMass)) * cellDensity * (units.s * timeStepSec)

				for i, row in enumerate(xrange(0, 2 * len(enzymeIds), 2)):
					countAxis = axesList[row]
					fluxAxis = axesList[row + 1]
					plotFlux = (((COUNTS_UNITS / VOLUME_UNITS) * enzymeFluxes[:, i]) / coefficient).asNumber(units.mmol / units.g / units.h)
					countAxis.plot(time_ / 3600., enzymeCounts[:, i], color = "b")
					fluxAxis.plot(time_ / 3600., plotFlux, color = "b")
				axesList[-2].plot(time_ / 3600., reactantCounts, color = "b")
				axesList[-1].plot(time_ / 3600., metCounts, color = "b")

			ylabels = ["menE", "menB", "menI", "menA", "ubiE", "CPD-12115", "Menaquinone"]
			for i, axis in enumerate(axesList[::2]):
				axis.set_xlim([0, time_[-1] / 3600.])
				axis.set_ylabel("%s" % ylabels[i], rotation = 0)
				whitePadSparklineAxis(axis, False)
			for axis in axesList[1::2]:
				axis.set_xlim([0, time_[-1] / 3600.])
				whitePadSparklineAxis(axis)
			axesList[-1].set_ylabel(ylabels[-1], rotation = 0)

			for axis in axesList:
				for i in xrange(len(patchStart)):
					width = time[patchEnd[i]] / 3600. - time[patchStart[i]] / 3600.
					if width <= 0.1:
						continue

					height = axis.get_ylim()[1] - axis.get_ylim()[0]
					axis.add_patch(patches.Rectangle((time[patchStart[i]] / 3600., axis.get_ylim()[0]), width, height, alpha = 0.25, color = "gray", linewidth = 0.))

			plt.subplots_adjust(hspace = 0.5, right = 0.95, bottom = 0.05, left = 0.11, top = 0.95)
			exportFigure(plt, plotOutDir, plotOutFileName + "__downstreamFluxes", metadata)


if __name__ == "__main__":
	Plot().cli()
