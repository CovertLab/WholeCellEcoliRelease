"""
Plots the downstream effects of the subgenerational expression of pabA and pabB
genes, used in Figure 4 of the paper.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/25/2019
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

ENZYME_RNA_IDS = ["EG10682_RNA[c]", "EG10683_RNA[c]"]
ENZYME_MONOMER_IDS = ["PABASYN-COMPII-MONOMER[c]", "PABASYN-COMPI-MONOMER[c]"]
ENZYME_COMPLEX_IDS = ["PABASYN-CPLX[c]", "PABSYNMULTI-CPLX[c]"]
ENZYME_REACTION_IDS = ["PABASYN-RXN__PABASYN-CPLX (reverse)", "PABASYN-RXN__PABSYNMULTI-CPLX (reverse)"]
METABOLITE_ID = "METHYLENE-THF[c]"

FLUX_LINEAR_THRESHOLD = 1e-3  # Flux values below this threshold will be plotted linearly

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
		rna_ids = sim_data.process.transcription.rnaData["id"]

		enzyme_rna_transcription_indexes = np.array([
			np.where(rna_ids == enzyme_rna_id)[0][0]
			for enzyme_rna_id in ENZYME_RNA_IDS
			])

		simOutDir = os.path.join(allDir[0], "simOut")
		bulk_molecules_reader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		fba_results_reader = TableReader(os.path.join(simOutDir, "FBAResults"))

		moleculeIDs = bulk_molecules_reader.readAttribute("objectNames")
		reactionIDs = np.array(fba_results_reader.readAttribute("reactionIDs"))

		enzyme_rna_count_indexes = np.array([
			moleculeIDs.index(enzyme_rna_id)
			for enzyme_rna_id in ENZYME_RNA_IDS
			])
		enzyme_monomer_indexes = np.array([
			moleculeIDs.index(enzyme_monomer_id)
			for enzyme_monomer_id in ENZYME_MONOMER_IDS
			])
		enzyme_complex_indexes = np.array([
			moleculeIDs.index(enzyme_complex_id)
			for enzyme_complex_id in ENZYME_COMPLEX_IDS
			])
		reaction_indexes = np.array([
			np.where(reactionIDs == reaction_id)[0][0]
			for reaction_id in ENZYME_REACTION_IDS
			])
		metabolite_index = moleculeIDs.index(METABOLITE_ID)

		# Initialize arrays
		time = []
		enzyme_rna_init_events = np.empty((0, len(ENZYME_RNA_IDS)))
		enzyme_rna_counts = np.empty((0, len(ENZYME_RNA_IDS)))
		enzyme_total_monomer_counts = np.empty((0, len(ENZYME_RNA_IDS)))
		enzyme_complex_counts = np.empty((0, len(ENZYME_COMPLEX_IDS)))
		enzyme_fluxes = np.empty((0, len(ENZYME_REACTION_IDS)))
		metabolite_counts = []

		cellMass = []
		dryMass = []
		timeStepSec = []
		generationTicks = []

		proteins_produced_per_gen = np.empty((0, len(ENZYME_RNA_IDS)))
		average_complex_counts_per_gen = []

		first_gen = True

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			mass_reader = TableReader(os.path.join(simOutDir, "Mass"))
			bulk_molecules_reader = TableReader(
				os.path.join(simOutDir, "BulkMolecules"))
			fba_results_reader = TableReader(os.path.join(simOutDir, "FBAResults"))
			rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

			time.extend(main_reader.readColumn("time").tolist())

			if first_gen:
				generationTicks.extend([time[0], time[-1]])
				first_gen = False
			else:
				generationTicks.append(time[-1])

			timeStepSec.extend(main_reader.readColumn("timeStepSec").tolist())
			cellMass.extend(mass_reader.readColumn("cellMass").tolist())
			dryMass.extend(mass_reader.readColumn("dryMass").tolist())

			rna_init_events_this_gen = rnap_data_reader.readColumn(
				"rnaInitEvent")[:, enzyme_rna_transcription_indexes]
			enzyme_rna_init_events = np.vstack((
				enzyme_rna_init_events, rna_init_events_this_gen))

			molecule_counts = bulk_molecules_reader.readColumn("counts")

			enzyme_rna_counts = np.vstack((
				enzyme_rna_counts,
				molecule_counts[:, enzyme_rna_count_indexes]))

			enzyme_monomer_counts_this_gen = molecule_counts[:, enzyme_monomer_indexes]
			enzyme_complex_counts_this_gen = molecule_counts[:, enzyme_complex_indexes]

			enzyme_total_monomer_counts_this_gen = (
					enzyme_monomer_counts_this_gen +
					enzyme_complex_counts_this_gen.sum(axis=1)[:, None])

			enzyme_total_monomer_counts = np.vstack((
				enzyme_total_monomer_counts,
				enzyme_total_monomer_counts_this_gen))
			enzyme_complex_counts = np.vstack((
				enzyme_complex_counts,
				enzyme_complex_counts_this_gen))
			proteins_produced_per_gen = np.vstack((
				proteins_produced_per_gen,
				(enzyme_total_monomer_counts_this_gen[-1, :] -
					enzyme_total_monomer_counts_this_gen[0, :])
				))
			average_complex_counts_per_gen.append(
				enzyme_complex_counts_this_gen.sum(axis=1).mean())

			metabolite_counts.extend(molecule_counts[:, metabolite_index])

			reactionFluxes = np.array(fba_results_reader.readColumn("reactionFluxes"))
			enzyme_fluxes = np.vstack((
				enzyme_fluxes, reactionFluxes[:, reaction_indexes]))

		# Sum reaction fluxes and convert units
		flux_conversion_coeff = (units.fg * np.array(dryMass)) / (units.fg * np.array(cellMass)) * (timeStepSec * units.s) * cellDensity
		enzyme_fluxes = (((COUNTS_UNITS / VOLUME_UNITS) * enzyme_fluxes.sum(axis=1)) / flux_conversion_coeff).asNumber(units.mmol / units.g / units.h)

		# Convert time to hours
		time = np.array(time)
		time_hours = time / 3600.

		# Add counts of complexed monomers to monomer counts
		enzyme_total_monomer_counts += enzyme_complex_counts.sum(axis=1)[:, None]

		# Plot
		plt.figure(figsize = (14, 8.5))
		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
		plt.suptitle(
			"4-amino-4-deoxychorismate synthase downstream effects",
			fontsize = FONTSIZE)
		pre_merge_colors = [color_cycle[0], color_cycle[2]]
		post_merge_color = color_cycle[3]

		# Define axes
		rna_init_axis = plt.subplot(6, 1, 1)
		rna_axis = plt.subplot(6, 1, 2, sharex=rna_init_axis)
		monomer_axis = plt.subplot(6, 1, 3, sharex=rna_init_axis)
		complex_axis = plt.subplot(6, 1, 4, sharex=rna_init_axis)
		flux_axis = plt.subplot(6, 1, 5, sharex=rna_init_axis)
		met_axis = plt.subplot(6, 1, 6, sharex=rna_init_axis)

		# Plot transcription initiation events
		rna_init_axis.set_prop_cycle(color=pre_merge_colors)
		rna_init_axis.plot(time_hours, enzyme_rna_init_events)
		rna_init_axis.set_ylabel("Transcription\nevents", fontsize = FONTSIZE, rotation = 0)
		rna_init_axis.yaxis.set_label_coords(-.12, 0.25)
		rna_init_axis.set_xlim([time_hours[0], time_hours[-1]])
		rna_init_axis.set_ylim([0, 1])
		whitePadSparklineAxis(rna_init_axis, xAxis = False)

		# Print average transcription frequency of each gene
		for rna_id, prob in zip(ENZYME_RNA_IDS,
				enzyme_rna_init_events.sum(axis=0)/len(proteins_produced_per_gen)):
			print("%s transcription frequency: %.3f"%(rna_id, prob))

		rna_axis.set_prop_cycle(color=pre_merge_colors)
		rna_axis.plot(time_hours, enzyme_rna_counts)
		rna_axis.set_ylabel("mRNA\ncounts", fontsize = FONTSIZE, rotation = 0)
		rna_axis.yaxis.set_label_coords(-.12, 0.25)
		rna_axis.set_ylim([0, np.max(enzyme_rna_counts)])
		whitePadSparklineAxis(rna_axis, xAxis = False)

		monomer_axis.set_prop_cycle(color=pre_merge_colors)
		monomer_axis.plot(time_hours, enzyme_total_monomer_counts)
		monomer_axis.set_ylabel("Protein monomer\ncounts", fontsize = FONTSIZE, rotation = 0)
		monomer_axis.yaxis.set_label_coords(-.12, 0.25)
		monomer_axis.set_ylim([0, np.max(enzyme_total_monomer_counts)])
		whitePadSparklineAxis(monomer_axis, xAxis = False)

		# Print average number of protein produced per generation
		for rna_id, count in zip(ENZYME_RNA_IDS,
				proteins_produced_per_gen.mean(axis=0)):
			print("%s average proteins produced per gen: %.2f" % (rna_id, count))

		complex_axis.plot(time_hours, enzyme_complex_counts.sum(axis=1), color=post_merge_color)
		complex_axis.set_ylabel("Protein complex\ncounts", fontsize = FONTSIZE, rotation = 0)
		complex_axis.yaxis.set_label_coords(-.12, 0.25)
		complex_axis.set_ylim([0, np.max(enzyme_complex_counts.sum(axis=1))])
		whitePadSparklineAxis(complex_axis, xAxis = False)

		# Print mean and std of average complex counts in each gen
		print("Complex counts average: %.2f" % (np.array(average_complex_counts_per_gen).mean(),))
		print("Complex counts std: %.2f" % (np.array(average_complex_counts_per_gen).std(),))

		flux_axis.plot(time_hours, enzyme_fluxes, color=post_merge_color)
		flux_axis.set_yscale("symlog", linthreshy=FLUX_LINEAR_THRESHOLD)
		flux_axis.set_ylabel("PABASYN-RXN\n(reverse)\ntotal flux\n(mmol/gDCW/hour)", fontsize = FONTSIZE, rotation = 0)
		flux_axis.yaxis.set_label_coords(-.12, 0.25)
		flux_axis.set_ylim([0, np.max(enzyme_fluxes)])
		whitePadSparklineAxis(flux_axis, xAxis=False)
		flux_axis.get_yaxis().set_tick_params(which='minor', size=0)
		flux_axis.get_xaxis().set_tick_params(which='minor', width=0)
		flux_max = flux_axis.get_ylim()[1]
		flux_axis.set_yticks([0, FLUX_LINEAR_THRESHOLD, flux_max])
		flux_axis.set_yticklabels(["0", "%0.0e"%(FLUX_LINEAR_THRESHOLD, ), "%.2f"%(flux_max, )])

		met_axis.plot(time_hours, metabolite_counts, color=post_merge_color)
		met_axis.set_ylabel("End product\ncounts", fontsize = FONTSIZE, rotation = 0)
		met_axis.yaxis.set_label_coords(-.12, 0.25)
		met_axis.set_xlabel("Time (hour)\ntickmarks at each new generation", fontsize = FONTSIZE)
		met_axis.set_ylim([0, np.max(metabolite_counts)])
		met_axis.set_xlim([time_hours[0], time_hours[-1]])
		whitePadSparklineAxis(met_axis)
		met_axis.set_yticklabels([0, "%0.1e" % met_axis.get_ylim()[1]])
		met_axis.set_xticks(np.array(generationTicks) / 3600.)
		xticklabels = np.repeat("     ", len(generationTicks))
		xticklabels[0] = "%0.2f" % (time_hours[0])
		xticklabels[-1] = "%0.2f" % (time_hours[-1])
		met_axis.set_xticklabels(xticklabels)

		# Add patches to indicate absence of complexes
		noComplexIndexes = np.where(np.array(enzyme_complex_counts.sum(axis=1)) == 0)[0]
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

		axesList = [rna_init_axis, rna_axis, monomer_axis, complex_axis, flux_axis, met_axis]
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
		met_axis.set_xticklabels([])
		met_axis.set_xlabel("")
		exportFigure(plt, plotOutDir, plotOutFileName + "__clean", "")
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
