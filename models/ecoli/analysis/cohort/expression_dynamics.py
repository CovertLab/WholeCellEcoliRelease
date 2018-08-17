from __future__ import absolute_import

import os
from itertools import izip
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.rdp import rdp
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

FROM_CACHE = False

GENS = np.arange(3, 9)


def mm2inch(value):
	return value * 0.0393701

def align_yaxis(ax1, v1, ax2, v2):
	"""adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
	_, y1 = ax1.transData.transform((0, v1))
	_, y2 = ax2.transData.transform((0, v2))
	inv = ax2.transData.inverted()
	_, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
	miny, maxy = ax2.get_ylim()
	ax2.set_ylim(miny+dy, maxy+dy)

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Check if basal sim
		sim_data = cPickle.load(open(simDataFile, "rb"))
		if sim_data.condition != "basal":
			print "Skipping - plot only runs for basal sim."
			return

		# Get all cells
		ap = AnalysisPaths(seedOutDir, cohort_plot = True)
		allDir = ap.get_cells(seed=[0], generation = GENS)
		n_gens = GENS.size
		if len(allDir) < n_gens:
			print "Skipping - particular seed and/or gens were not simulated."
			return

		# Get all ids reqiured
		ids_complexation = sim_data.process.complexation.moleculeNames # Complexes of proteins, and protein monomers
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes # Only complexes
		ids_equilibrium = sim_data.process.equilibrium.moleculeNames # Complexes of proteins + small molecules, small molecules, protein monomers
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes # Only complexes
		ids_translation = sim_data.process.translation.monomerData["id"].tolist() # Only protein monomers
		ids_transcription = sim_data.process.transcription.rnaData["id"].tolist()

		data_50s = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)
		data_30s = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)
		ribosome_subunit_ids = data_50s["subunitIds"].tolist() + data_30s["subunitIds"].tolist()
		ribosome_subunit_stoich = np.hstack((data_50s["subunitStoich"],data_30s["subunitStoich"]))
		data_rnap = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.rnapFull)
		rnap_subunit_ids = data_rnap["subunitIds"].tolist()
		rnap_subunit_stoich = data_rnap["subunitStoich"]

		# Stoich matrices
		complexStoich = sim_data.process.complexation.stoichMatrixMonomers()
		equilibriumStoich = sim_data.process.equilibrium.stoichMatrixMonomers()

		first_build = True

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomerData['id'].size
		n_sims = ap.n_generation

		ratioFinalToInitialCountMultigen = np.zeros((n_gens, n_monomers), dtype = np.float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)

		# protein_index_of_interest_full = np.zeros((n_gens, n_monomers), dtype = np.bool)

		if not FROM_CACHE:
			print "Re-running - not using cache"
			for gen_idx, simDir in enumerate(allDir):
				simOutDir = os.path.join(simDir, "simOut")

				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				## READ DATA ##
				# Read in bulk ids and counts
				bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

				if first_build:
					print "Running first build code"
					moleculeIds = bulkMolecules.readAttribute("objectNames")

					complexationIdx = np.array([moleculeIds.index(x) for x in ids_complexation]) # Complexes of proteins, and protein monomers
					complexation_complexesIdx = np.array([moleculeIds.index(x) for x in ids_complexation_complexes]) # Only complexes
					equilibriumIdx = np.array([moleculeIds.index(x) for x in ids_equilibrium]) # Complexes of proteins + small molecules, small molecules, protein monomers
					equilibrium_complexesIdx = np.array([moleculeIds.index(x) for x in ids_equilibrium_complexes]) # Only complexes
					translationIdx = np.array([moleculeIds.index(x) for x in ids_translation]) # Only protein monomers
					transcriptionIdx = np.array([moleculeIds.index(x) for x in ids_transcription]) # Only protein rnas
					ribosomeIdx = np.array([moleculeIds.index(x) for x in ribosome_subunit_ids])
					rnapIdx = np.array([moleculeIds.index(x) for x in rnap_subunit_ids])

					cPickle.dump(complexationIdx, open(os.path.join(plotOutDir,"complexationIdx.pickle"), "wb"))
					cPickle.dump(complexation_complexesIdx, open(os.path.join(plotOutDir,"complexation_complexesIdx.pickle"), "wb"))
					cPickle.dump(equilibriumIdx, open(os.path.join(plotOutDir,"equilibriumIdx.pickle"), "wb"))
					cPickle.dump(equilibrium_complexesIdx, open(os.path.join(plotOutDir,"equilibrium_complexesIdx.pickle"), "wb"))
					cPickle.dump(translationIdx, open(os.path.join(plotOutDir,"translationIdx.pickle"), "wb"))
					cPickle.dump(transcriptionIdx, open(os.path.join(plotOutDir,"transcriptionIdx.pickle"), "wb"))
					cPickle.dump(ribosomeIdx, open(os.path.join(plotOutDir,"ribosomeIdx.pickle"), "wb"))
					cPickle.dump(rnapIdx, open(os.path.join(plotOutDir,"rnapIdx.pickle"), "wb"))

					first_build = False

				bulkCounts = bulkMolecules.readColumn("counts")
				bulkMolecules.close()

				# Dissociate protein-protein complexes
				bulkCounts[:, complexationIdx] += np.dot(complexStoich, bulkCounts[:, complexation_complexesIdx].transpose() * -1).transpose().astype(np.int)

				# Dissociate protein-small molecule complexes
				bulkCounts[:, equilibriumIdx] += np.dot(equilibriumStoich, bulkCounts[:, equilibrium_complexesIdx].transpose() * -1).transpose().astype(np.int)

				# Load unique molecule data for RNAP and ribosomes
				uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
				rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
				nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
				nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
				uniqueMoleculeCounts.close()

				# Add subunits from RNAP and ribosomes
				ribosomeSubunitCounts = (nActiveRibosome.reshape((nActiveRibosome.size,1)) * ribosome_subunit_stoich.reshape((1,ribosome_subunit_stoich.size)))
				rnapSubunitCounts = (nActiveRnaPoly.reshape((nActiveRnaPoly.size,1)) * rnap_subunit_stoich.reshape((1,rnap_subunit_stoich.size)))

				bulkCounts[:, ribosomeIdx] += ribosomeSubunitCounts.astype(np.int64)
				bulkCounts[:, rnapIdx] += rnapSubunitCounts.astype(np.int64)

				# Get protein monomer counts for calculations now that all complexes are dissociated
				proteinMonomerCounts = bulkCounts[:, translationIdx]
				ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(np.float) + 1)
				ratioFinalToInitialCountMultigen[gen_idx, :] = ratioFinalToInitialCount

			# cPickle.dump(protein_index_of_interest_full, open(os.path.join(plotOutDir,"protein_index_of_interest_full.pickle"), "wb"))
			cPickle.dump(ratioFinalToInitialCountMultigen, open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "wb"))

		# protein_index_of_interest_full = cPickle.load(open(os.path.join(plotOutDir,"protein_index_of_interest_full.pickle"), "rb"))
		# protein_index_of_interest = np.where(protein_index_of_interest_full.all(axis = 0))[0]

		ratioFinalToInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "rb"))
		complexationIdx = cPickle.load(open(os.path.join(plotOutDir,"complexationIdx.pickle"), "rb"))
		complexation_complexesIdx = cPickle.load(open(os.path.join(plotOutDir,"complexation_complexesIdx.pickle"), "rb"))
		equilibriumIdx = cPickle.load(open(os.path.join(plotOutDir,"equilibriumIdx.pickle"), "rb"))
		equilibrium_complexesIdx = cPickle.load(open(os.path.join(plotOutDir,"equilibrium_complexesIdx.pickle"), "rb"))
		translationIdx = cPickle.load(open(os.path.join(plotOutDir,"translationIdx.pickle"), "rb"))
		transcriptionIdx = cPickle.load(open(os.path.join(plotOutDir,"transcriptionIdx.pickle"), "rb"))
		ribosomeIdx = cPickle.load(open(os.path.join(plotOutDir,"ribosomeIdx.pickle"), "rb"))
		rnapIdx = cPickle.load(open(os.path.join(plotOutDir,"rnapIdx.pickle"), "rb"))

		protein_index_of_interest = np.where(np.logical_and(ratioFinalToInitialCountMultigen > 1.6, ratioFinalToInitialCountMultigen < 2.4).all(axis = 0))[0]
		first_gen_flat = ratioFinalToInitialCountMultigen[0,:] < 1.1
		second_gen_burst = ratioFinalToInitialCountMultigen[1,:] > 10
		rest_of_gens_decline = (ratioFinalToInitialCountMultigen[2:,:] < 1.1).all(axis=0)
		logic_filter = np.logical_and.reduce((first_gen_flat, second_gen_burst, rest_of_gens_decline))
		protein_index_of_interest_burst = np.where(logic_filter)[0]
		try: # try block expects particular proteins to plot
			protein_index_of_interest = protein_index_of_interest[:5]
			protein_idx = protein_index_of_interest[4]
			protein_idx_burst = protein_index_of_interest_burst[2]
		except Exception as exc:
			print "Error: %s" % exc
			return

		fig, axesList = plt.subplots(ncols = 2, nrows = 2, sharex = True)
		expProtein_axis = axesList[0,0]
		expRna_axis = axesList[1,0]
		burstProtein_axis = axesList[0,1]
		burstRna_axis = axesList[1,1]

		# expProteinFold_axis = expProtein_axis.twinx()
		# expProteinFold_axis.spines["bottom"].set_visible(False)
		# expProteinFold_axis.spines["top"].set_visible(False)
		# expProteinFold_axis.spines["left"].set_visible(False)
		# expProteinFold_axis.tick_params(bottom = "off")
		# expProteinFold_axis.tick_params(axis = "x", labelbottom='off')
		# expProteinFold_axis.spines["right"].set_position(("outward", 10))
		# expProteinFold_axis.tick_params(which = "both", direction = "out")
		# expProteinFold_axis.set_ylabel("Fold change", fontsize=9)

		# burstProteinFold_axis = burstProtein_axis.twinx()
		# burstProteinFold_axis.spines["bottom"].set_visible(False)
		# burstProteinFold_axis.spines["top"].set_visible(False)
		# burstProteinFold_axis.spines["left"].set_visible(False)
		# burstProteinFold_axis.tick_params(bottom = "off")
		# burstProteinFold_axis.tick_params(axis = "x", labelbottom='off')
		# burstProteinFold_axis.spines["right"].set_position(("outward", 10))
		# burstProteinFold_axis.tick_params(which = "both", direction = "out")
		# burstProteinFold_axis.set_ylabel("Fold change", fontsize=9)

		mult = 3
		fig_width = mm2inch(80) * mult
		fig_height = mm2inch(50) * mult

		fig.set_figwidth(fig_width)
		fig.set_figheight(fig_height)
		firstLine = True
		firstLineInit = None
		firstLineInitRna = None
		firstLineInit_burst = None
		firstLineInitRna_burst = None

		time_eachGen = []
		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			time_eachGen.append(time[0])
			if gen_idx == 0:
				startTime = time[0]

			## READ DATA ##
			# Read in bulk ids and counts
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkCounts = bulkMolecules.readColumn("counts")
			bulkMolecules.close() # NOTE (John): .close() doesn't currently do anything

			# Dissociate protein-protein complexes
			bulkCounts[:, complexationIdx] += np.dot(complexStoich, bulkCounts[:, complexation_complexesIdx].transpose() * -1).transpose().astype(np.int)

			# Dissociate protein-small molecule complexes
			bulkCounts[:, equilibriumIdx] += np.dot(equilibriumStoich, bulkCounts[:, equilibrium_complexesIdx].transpose() * -1).transpose().astype(np.int)

			# Load unique molecule data for RNAP and ribosomes
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
			nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
			uniqueMoleculeCounts.close()

			# Add subunits from RNAP and ribosomes
			ribosomeSubunitCounts = (nActiveRibosome.reshape((nActiveRibosome.size,1)) * ribosome_subunit_stoich.reshape((1,ribosome_subunit_stoich.size)))
			rnapSubunitCounts = (nActiveRnaPoly.reshape((nActiveRnaPoly.size,1)) * rnap_subunit_stoich.reshape((1,rnap_subunit_stoich.size)))

			bulkCounts[:, ribosomeIdx] += ribosomeSubunitCounts.astype(np.int64)
			bulkCounts[:, rnapIdx] += rnapSubunitCounts.astype(np.int64)

			# Get protein monomer counts for calculations now that all complexes are dissociated
			proteinMonomerCounts = bulkCounts[:, translationIdx]
			rnaMonomerCounts = bulkCounts[:, transcriptionIdx]

			if firstLine:
				firstLineInit = float(proteinMonomerCounts[:, protein_idx][0])
				firstLineInitRna = float(rnaMonomerCounts[:, sim_data.relation.rnaIndexToMonomerMapping][:,protein_idx][0])
				firstLineInit_burst = float(proteinMonomerCounts[:, protein_idx_burst][0])
				firstLineInitRna_burst = float(rnaMonomerCounts[:, sim_data.relation.rnaIndexToMonomerMapping][:,protein_idx_burst][0])
				firstLine = False

			LINEWIDTH = 1

			time_minutes = time / 60.

			EXP_COLOR = 'red'
			BURST_COLOR = 'blue'

			axes = (
				expProtein_axis,
				burstProtein_axis,
				expRna_axis,
				burstRna_axis
				)
			counts = (
				proteinMonomerCounts[:, protein_idx],
				proteinMonomerCounts[:, protein_idx_burst],
				rnaMonomerCounts[:, sim_data.relation.rnaIndexToMonomerMapping][:,protein_idx],
				rnaMonomerCounts[:, sim_data.relation.rnaIndexToMonomerMapping][:,protein_idx_burst]
				)
			line_color = (
				EXP_COLOR,
				BURST_COLOR,
				EXP_COLOR,
				BURST_COLOR,
				)
			count_min = ( # better to acquire programatically, but would require loading data twice
				600,
				0,
				0,
				0
				)
			count_scale = ( # same as above
				2200 - 600,
				30 - 0,
				7 - 0,
				1 - 0
				)

			# These are *approximate* estimates of the axes sizes, using the
			# size of the figure plus the fact that the subplots are 2x2.
			# This is good enough since we primarily need the aspect ratio;
			# however there is probably a programmatic way to get this info
			# from the axes objects themselves.
			axes_width = fig_width / 2
			axes_height = fig_height / 2

			# No easy way to know how long the total set of simulations
			# will be without rewriting a lot of code, so assume that
			# the total time is roughly the time of the current generation
			# multiplied by the total number of generations.
			rescaled_time = (time_minutes - time_minutes.min())/(
				(time_minutes.max() - time_minutes.min()) * n_gens
				)

			for (ax, c, lc, cm, cs) in izip(axes, counts, line_color, count_min, count_scale):
				rescaled_counts = (c.astype(np.float64) - cm)/cs

				# Roughly rescale the data into the plotted dimensions for
				# better point downsampling
				points = np.column_stack([
					rescaled_time * axes_width,
					rescaled_counts * axes_height
					])

				RDP_THRESHOLD = 1e-5

				keep = rdp(points, RDP_THRESHOLD)

				x = time_minutes[keep]
				y = c[keep]

				ax.plot(x, y, color = lc, linewidth = LINEWIDTH)

			# expProteinFold_axis.plot(time_minutes, proteinMonomerCounts[:, protein_idx] / firstLineInit, alpha = 0.,color = "red")
			# burstProteinFold_axis.plot(time_minutes, proteinMonomerCounts[:, protein_idx_burst] / firstLineInit_burst, alpha = 0., color="red")

		expProtein_axis.set_title("Exponential dynamics: {}".format(sim_data.process.translation.monomerData['id'][protein_idx][:-3]), fontsize=9)
		burstProtein_axis.set_title("Sub-generational dynamics: {}".format(sim_data.process.translation.monomerData['id'][protein_idx_burst][:-3]), fontsize=9)
		expProtein_axis.set_ylabel("Protein\ncount", rotation=0, fontsize=9)
		expRna_axis.set_ylabel("mRNA\ncount", rotation=0, fontsize=9)

		# align_yaxis(expProtein_axis, firstLineInit, expProteinFold_axis, 1)
		# expProteinFold_axis.set_yticks([expProteinFold_axis.get_ylim()[0], 1., expProteinFold_axis.get_ylim()[1]])
		# expProteinFold_axis.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

		# align_yaxis(burstProtein_axis, firstLineInit_burst, burstProteinFold_axis, 1)
		# burstProteinFold_axis.set_yticks([burstProteinFold_axis.get_ylim()[0], 1., burstProteinFold_axis.get_ylim()[1]])

		time_eachGen.append(time[-1])
		time_eachGen = np.array(time_eachGen)

		expProtein_axis.set_xlim([startTime / 60., time[-1] / 60.])
		burstProtein_axis.set_xlim([startTime / 60., time[-1] / 60.])
		expRna_axis.set_xlim([startTime / 60., time[-1] / 60.])
		burstRna_axis.set_xlim([startTime / 60., time[-1] / 60.])

		whitePadSparklineAxis(expProtein_axis, False)
		whitePadSparklineAxis(burstProtein_axis, False)
		whitePadSparklineAxis(expRna_axis)
		whitePadSparklineAxis(burstRna_axis)

		expRna_axis.set_xticks(time_eachGen / 60.)
		burstRna_axis.set_xticks(time_eachGen / 60.)
		xlabel = GENS.tolist()
		xlabel.append(GENS[-1] + 1)
		expRna_axis.set_xticklabels(xlabel)
		burstRna_axis.set_xticklabels(xlabel)

		burstRna_axis.set_xlabel("Time (gens)", fontsize = 9)
		expRna_axis.set_xlabel("Time (gens)", fontsize = 9)

		axesList = axesList.flatten().tolist()
		# axesList.append(expProteinFold_axis)
		# axesList.append(burstProteinFold_axis)
		for axes in axesList:
			for tick in axes.xaxis.get_major_ticks():
				tick.label.set_fontsize(9)
			for tick in axes.yaxis.get_major_ticks():
				tick.label.set_fontsize(9)

		plt.subplots_adjust(wspace = 0.2, hspace = 0.15)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		for axes in axesList:
			axes.set_xlabel("")
			axes.set_ylabel("")
			axes.set_title("")
			axes.set_xticklabels([])
			axes.set_yticklabels([])

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
