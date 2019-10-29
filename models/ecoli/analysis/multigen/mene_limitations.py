"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/12/2017
"""

from __future__ import absolute_import

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
from models.ecoli.analysis import multigenAnalysisPlot

FONTSIZE = 6
LABELSIZE = 6

ENZYME_COMPLEX_ID = "MENE-CPLX[c]"
ENZYME_MONOMER_ID = "O-SUCCINYLBENZOATE-COA-LIG-MONOMER[c]"
ENZYME_RNA_ID = "EG12437_RNA[c]"
ENZYME_REACTION_ID = "O-SUCCINYLBENZOATE-COA-LIG-RXN"
METABOLITE_IDS = ["REDUCED-MENAQUINONE[c]", "CPD-12115[c]"]

def clearLabels(axis):
	axis.set_yticklabels([])
	axis.set_ylabel("")

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		sim_data = cPickle.load(open(simDataFile, "rb"))
		cellDensity = sim_data.constants.cellDensity
		rnaIds = sim_data.process.transcription.rnaData["id"]

		enzyme_rna_transcription_index = np.where(rnaIds == ENZYME_RNA_ID)[0][0]

		simOutDir = os.path.join(allDir[0], "simOut")
		bulk_molecules_reader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		fba_results_reader = TableReader(os.path.join(simOutDir, "FBAResults"))

		moleculeIDs = bulk_molecules_reader.readAttribute("objectNames")
		reactionIDs = np.array(fba_results_reader.readAttribute("reactionIDs"))

		enzyme_complex_index = moleculeIDs.index(ENZYME_COMPLEX_ID)
		enzyme_monomer_index = moleculeIDs.index(ENZYME_MONOMER_ID)
		enzyme_rna_counts_index = moleculeIDs.index(ENZYME_RNA_ID)
		metabolite_indexes = [moleculeIDs.index(x) for x in METABOLITE_IDS]
		reaction_index = np.where(reactionIDs == ENZYME_REACTION_ID)[0][0]

		# Initialize arrays
		time = []
		enzyme_fluxes = []
		enzyme_complex_counts = []
		enzyme_monomer_counts = []
		enzyme_rna_counts = []
		enzyme_rna_init_events = []
		metabolite_counts = np.empty((0, len(METABOLITE_IDS)))

		cellMass = []
		dryMass = []
		timeStepSec = []
		generationTicks = [0.]

		n_transcription_init_events_per_gen = []
		enzyme_complex_avg_counts = []

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			mass_reader = TableReader(os.path.join(simOutDir, "Mass"))
			bulk_molecules_reader = TableReader(
				os.path.join(simOutDir, "BulkMolecules"))
			fba_results_reader = TableReader(os.path.join(simOutDir, "FBAResults"))
			rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

			time.extend(main_reader.readColumn("time").tolist())
			generationTicks.append(time[-1])

			timeStepSec.extend(main_reader.readColumn("timeStepSec").tolist())
			cellMass.extend(mass_reader.readColumn("cellMass").tolist())
			dryMass.extend(mass_reader.readColumn("dryMass").tolist())

			molecule_counts = bulk_molecules_reader.readColumn("counts")
			enzyme_monomer_counts.extend(molecule_counts[:, enzyme_monomer_index].tolist())
			enzyme_rna_counts.extend(molecule_counts[:, enzyme_rna_counts_index].tolist())

			enzyme_complex_counts_this_gen = molecule_counts[:, enzyme_complex_index]
			enzyme_complex_counts.extend(enzyme_complex_counts_this_gen.tolist())
			enzyme_complex_avg_counts.append(np.mean(enzyme_complex_counts_this_gen))

			metabolite_counts = np.vstack((
				metabolite_counts, molecule_counts[:, metabolite_indexes]))

			reactionFluxes = np.array(fba_results_reader.readColumn("reactionFluxes"))
			enzyme_fluxes.extend(reactionFluxes[:, reaction_index].tolist())

			rna_init_events_this_gen = rnap_data_reader.readColumn("rnaInitEvent")[:, enzyme_rna_transcription_index]

			enzyme_rna_init_events.extend(rna_init_events_this_gen.tolist())
			n_transcription_init_events_per_gen.append(np.sum(rna_init_events_this_gen))

		coefficient = (units.fg * np.array(dryMass)) / (units.fg * np.array(cellMass)) * (timeStepSec * units.s) * cellDensity
		enzyme_fluxes = (((COUNTS_UNITS / VOLUME_UNITS) * enzyme_fluxes) / coefficient).asNumber(units.mmol / units.g / units.h)

		# Convert time to hours
		time = np.array(time)
		time_hours = time / 3600.

		# Plot
		plt.figure(figsize = (11, 8.5))
		plt.suptitle("O-succinylbenzoate-CoA ligase downstream behaviors", fontsize = FONTSIZE)

		# Define axes
		rnaInitAxis = plt.subplot(6, 1, 1)
		rnaAxis = plt.subplot(6, 1, 2, sharex = rnaInitAxis)
		monomerAxis = plt.subplot(6, 1, 3, sharex = rnaInitAxis)
		complexAxis = plt.subplot(6, 1, 4, sharex = rnaInitAxis)
		fluxAxis = plt.subplot(6, 1, 5, sharex = rnaInitAxis)
		metAxis = plt.subplot(6, 1, 6)

		# Plot transcription initiation events
		rnaInitAxis.plot(time_hours, enzyme_rna_init_events, c = "b")
		rnaInitAxis.set_ylabel(r"$menE$" + "\n transcription\nevents", fontsize = FONTSIZE, rotation = 0)
		rnaInitAxis.yaxis.set_label_coords(-.1, 0.25)
		rnaInitAxis.set_xlim([time_hours[0], time_hours[-1]])
		whitePadSparklineAxis(rnaInitAxis, xAxis = False)
		rnaInitAxis.set_yticks([0, 1])

		rnaAxis.plot(time_hours, enzyme_rna_counts, c = "b")
		rnaAxis.set_ylabel("menE mRNA\ncounts", fontsize = FONTSIZE, rotation = 0)
		rnaAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(rnaAxis, xAxis = False)
		rnaAxis.set_yticks([0, max(enzyme_rna_counts)])

		monomerAxis.plot(time_hours, enzyme_monomer_counts, c = "b")
		monomerAxis.set_ylabel("MenE monomer\ncounts", fontsize = FONTSIZE, rotation = 0)
		monomerAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(monomerAxis, xAxis = False)
		monomerAxis.set_yticks([0, 4, max(enzyme_monomer_counts)])

		complexAxis.plot(time_hours, enzyme_complex_counts, c = "b")
		complexAxis.set_ylabel("MenE tetramer\ncounts", fontsize = FONTSIZE, rotation = 0)
		complexAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(complexAxis, xAxis = False)
		complexAxis.set_yticks([0, max(enzyme_complex_counts)])

		fluxAxis.plot(time_hours, enzyme_fluxes, c = "b")
		fluxAxis.set_ylabel("SUCBZL flux\n(mmol/gDCW/hour)", fontsize = FONTSIZE, rotation = 0)
		fluxAxis.yaxis.set_label_coords(-.1, 0.25)
		whitePadSparklineAxis(fluxAxis, xAxis = False)
		fluxAxis.set_yticks([min(enzyme_fluxes), max(enzyme_fluxes)])

		metAxis.plot(time_hours, np.sum(metabolite_counts, axis = 1), c = "b")
		metAxis.set_ylabel("End product\ncounts", fontsize = FONTSIZE, rotation = 0)
		metAxis.yaxis.set_label_coords(-.1, 0.25)
		metAxis.set_xlabel("Time (hour)\ntickmarks at each new generation", fontsize = FONTSIZE)
		metAxis.set_ylim([metAxis.get_ylim()[0] * 0.2, metAxis.get_ylim()[1]])
		metAxis.set_xlim([time_hours[0], time_hours[-1]])
		whitePadSparklineAxis(metAxis)
		metAxis.set_yticklabels(["%0.1e" % metAxis.get_ylim()[0], "%0.1e" % metAxis.get_ylim()[1]])
		metAxis.set_xticks(np.array(generationTicks) / 3600.)
		xticklabels = np.repeat("     ", len(generationTicks))
		xticklabels[0] = "0"
		xticklabels[-1] = "%0.2f" % (time_hours[-1])
		metAxis.set_xticklabels(xticklabels)

		noComplexIndexes = np.where(np.array(enzyme_complex_counts) == 0)[0]
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
				width = time_hours[patchEnd[i]] - time_hours[patchStart[i]]
				if width <= 0.1:
					continue

				height = axis.get_ylim()[1] - axis.get_ylim()[0]
				axis.add_patch(patches.Rectangle((time_hours[patchStart[i]], axis.get_ylim()[0]), width, height, alpha = 0.25, color = "gray", linewidth = 0.))

		plt.subplots_adjust(hspace = 0.5, right = 0.9, bottom = 0.1, left = 0.15, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Get clean version of plot
		for a in axesList:
			clearLabels(a)

		plt.suptitle("")
		metAxis.set_xticklabels([])
		metAxis.set_xlabel("")
		exportFigure(plt, plotOutDir, plotOutFileName + "__clean", "")
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
