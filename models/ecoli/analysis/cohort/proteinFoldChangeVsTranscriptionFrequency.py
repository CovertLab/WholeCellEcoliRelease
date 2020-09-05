from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from six.moves import cPickle

from wholecell.utils import units

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

CLOSE_TO_DOUBLE = 0.1
FONT_SIZE = 9


def mm2inch(value):
	return value * 0.0393701

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		# Get all ids reqiured
		sim_data = cPickle.load(open(simDataFile, "rb"))
		ids_complexation = sim_data.process.complexation.molecule_names # Complexe of proteins, and protein monomers
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes # Only complexes
		ids_equilibrium = sim_data.process.equilibrium.molecule_names # Complexes of proteins + small molecules, small molecules, protein monomers
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes # Only complexes
		ids_translation = sim_data.process.translation.monomer_data["id"].tolist() # Only protein monomers

		# ids_ribosome =
		data_50s = sim_data.process.complexation.get_monomers(sim_data.molecule_ids.s50_full_complex)
		data_30s = sim_data.process.complexation.get_monomers(sim_data.molecule_ids.s30_full_complex)
		ribosome_subunit_ids = data_50s["subunitIds"].tolist() + data_30s["subunitIds"].tolist()
		ribosome_subunit_stoich = np.hstack((data_50s["subunitStoich"],data_30s["subunitStoich"]))

		data_rnap = sim_data.process.complexation.get_monomers(sim_data.molecule_ids.full_RNAP)
		rnap_subunit_ids = data_rnap["subunitIds"].tolist()
		rnap_subunit_stoich = data_rnap["subunitStoich"]

		# Get all cells
		ap = AnalysisPaths(seedOutDir, cohort_plot = True)
		allDir = ap.get_cells(seed = [0])

		first_build = True

		# Pre-allocate variables. Rows = Generations, Cols = Monomers
		n_monomers = sim_data.process.translation.monomer_data['id'].size
		n_sims = ap.n_generation

		monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = np.bool)
		ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = np.float)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
		monomerCountInitialMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
		cellMassInitialMultigen = np.zeros(n_sims, dtype = np.float)

		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			## READ DATA ##
			# Read in bulk ids and counts
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

			if first_build:
				moleculeIds = bulkMolecules.readAttribute("objectNames")

				complexationIdx = np.array([moleculeIds.index(x) for x in ids_complexation]) # Complexe of proteins, and protein monomers
				complexation_complexesIdx = np.array([moleculeIds.index(x) for x in ids_complexation_complexes]) # Only complexes
				equilibriumIdx = np.array([moleculeIds.index(x) for x in ids_equilibrium]) # Complexes of proteins + small molecules, small molecules, protein monomers
				equilibrium_complexesIdx = np.array([moleculeIds.index(x) for x in ids_equilibrium_complexes]) # Only complexes
				translationIdx = np.array([moleculeIds.index(x) for x in ids_translation]) # Only protein monomers

				ribosomeIdx = np.array([moleculeIds.index(x) for x in ribosome_subunit_ids])
				rnapIdx = np.array([moleculeIds.index(x) for x in rnap_subunit_ids])

				first_build = False

			bulkCounts = bulkMolecules.readColumn("counts")
			bulkMolecules.close()

			# Dissociate protein-protein complexes
			bulkCounts[:, complexationIdx] += np.dot(sim_data.process.complexation.stoich_matrix_monomers(), bulkCounts[:, complexation_complexesIdx].transpose() * -1).transpose()

			# Dissociate protein-small molecule complexes
			bulkCounts[:, equilibriumIdx] += np.dot(sim_data.process.equilibrium.stoich_matrix_monomers(), bulkCounts[:, equilibrium_complexesIdx].transpose() * -1).transpose()

			# Load unique molecule data for RNAP and ribosomes
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
			rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_RNAP')
			nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
			uniqueMoleculeCounts.close()

			# Add subunits from RNAP and ribosomes
			ribosomeSubunitCounts = (nActiveRibosome.reshape((nActiveRibosome.size,1)) * ribosome_subunit_stoich.reshape((1,ribosome_subunit_stoich.size)))
			rnapSubunitCounts = (nActiveRnaPoly.reshape((nActiveRnaPoly.size,1)) * rnap_subunit_stoich.reshape((1,rnap_subunit_stoich.size)))

			bulkCounts[:, ribosomeIdx] += ribosomeSubunitCounts
			bulkCounts[:, rnapIdx] += rnapSubunitCounts

			# Get protein monomer counts for calculations now that all complexes are dissociated
			proteinMonomerCounts = bulkCounts[:, translationIdx]

			## CALCULATIONS ##
			# Calculate if monomer exists over course of cell cycle
			monomerExist = proteinMonomerCounts.sum(axis=0) > 1

			# Calculate if monomer comes close to doubling
			ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(np.float) + 1)
			# monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)

			# Load transcription initiation event data
			rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

			# Map transcription initiation events to monomers
			initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rna_index_to_monomer_mapping]

			# Load cell mass
			cellMassInitial = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")[0]

			# Log data
			monomerExistMultigen[gen_idx,:] = monomerExist
			ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer
			monomerCountInitialMultigen[gen_idx,:] = proteinMonomerCounts[0,:]
			cellMassInitialMultigen[gen_idx] = cellMassInitial

		cellMassInitialMultigen = units.fg * cellMassInitialMultigen

		existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
		averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen#.mean(axis=0)
		averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)

		averageInitiationEventsPerMonomer = np.tile(averageInitiationEventsPerMonomer, (6,1))


		mws = sim_data.getter.get_mass(sim_data.process.translation.monomer_data['id'])
		monomerInitialMasses = (mws * monomerCountInitialMultigen / sim_data.constants.n_Avogadro)

		# np.tile(cellMassInitialMultigen.asNumber().reshape((1,10)), (n_monomers,1))

		# initialMassFractions = monomerInitialMasses.asNumber(units.fg).transpose() / np.tile(cellMassInitialMultigen.asNumber().reshape((1,10)), (n_monomers,1))
		# averageInitialMassFractions = initialMassFractions.mean(axis = 1)
		# avgMonomerInitialMass = monomerInitialMasses.asNumber(units.fg).mean(axis=0)
		# avgMonomerInitialMassFraction = avgMonomerInitialMass / avgMonomerInitialMass.sum()

		# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		# probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		# probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		# for idx, burstSize in enumerate(uniqueBurstSizes):
		# 	mask = initiationEventsPerMonomerMultigen == burstSize
		# 	mask_sum = mask.sum()
		# 	probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
		# 	probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())


		# fig, axesList = plt.subplots(4,1)

		mult = 3
		fig = plt.figure()
		fig.set_figwidth(mm2inch(80) * mult)
		fig.set_figheight(mm2inch(50) * mult)

		scatterAxis = plt.subplot2grid((3,4), (0, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		# scatterAxis.axhline(1.0, linewidth=0.5, color='black', linestyle="--", xmin = 0.5, xmax = 1.)
		# scatterAxis.axhline(2.0, linewidth=0.5, color='black', linestyle="--", xmin = 0.5, xmax = 1.)
		# xhistAxis = plt.subplot2grid((4,5), (0,0), colspan=3, sharex = scatterAxis)
		yhistAxis = plt.subplot2grid((3,4), (0,3), rowspan=3)#, sharey = scatterAxis)
		# yhistAxis.axhline(1.0, linewidth=1.0, color='black', linestyle = 'dotted')
		# yhistAxis.axhline(2.0, linewidth=1.0, color='black', linestyle = 'dotted')
		#yhistAxis_2 = plt.subplot2grid((4,5), (1,4), rowspan=3, sharey = scatterAxis)
		#yhistAxis_2.axhline(1.0, linewidth=0.5, color='black', linestyle="--")
		#yhistAxis_2.axhline(2.0, linewidth=0.5, color='black', linestyle="--")

		# xhistAxis.xaxis.set_visible(False)
		#yhistAxis.yaxis.set_visible(False)


		smallBurst = averageInitiationEventsPerMonomer <= 1.
		# scatterAxis.set_xlim([1e-1, 1e3])
		# ----> scatterAxis.set_ylim([0.7, 100])

		# scatterAxis.semilogx(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "red", alpha = 0.9, lw = 0.)#, s = 5)
		# scatterAxis.semilogx(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "blue", alpha = 0.9, lw = 0.)#, s = 5)
		## scatterAxis.semilogx(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "green", alpha = 0.9, lw = 0.)#, s = 5)
		## scatterAxis.semilogx(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "blue", alpha = 0.9, lw = 0.)#, s = 5)

		scatterAxis.loglog(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "blue", alpha = 0.5, lw = 0.)#, s = 5)
		scatterAxis.loglog(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "red", alpha = 0.5, lw = 0.)#, s = 5)

		scatterAxis.set_ylabel("Fold change per protein\nin each generation ({} generations)".format(ap.n_generation), fontsize = FONT_SIZE)
		scatterAxis.set_xlabel("Average number of transcription events\nper protein per generation ({} generations)".format(ap.n_generation), fontsize = FONT_SIZE)

		# lims = yhistAxis.get_ylim()
		# step = (lims[1] - lims[0]) / 125
		# bins = np.arange(lims[0], lims[1] + step, step)

		# mass_in_binrange = np.zeros(bins.size-1, dtype=np.float)
		# for i in range(len(bins) - 1):
		# 	in_bin_range = np.logical_and(averageFoldChangePerMonomer > bins[i], averageFoldChangePerMonomer < bins[i+1])
		# 	mass_in_binrange[i] = avgMonomerInitialMassFraction[in_bin_range].sum()

		#yhistAxis_2.barh(bottom = bins[:-1], width = mass_in_binrange, height=(lims[1] - lims[0]) / 125, color = "white")
		#yhistAxis_2.set_xlim([0., 1.])
		#yhistAxis_2.yaxis.set_label_position("right")
		#yhistAxis_2.set_xlabel("Fraction of\nproteome mass", fontsize = FONT_SIZE)
		#scatterAxis.set_xlim([-10., 1000.])

		# yhistAxis.hist(averageFoldChangePerMonomer[~smallBurst], histtype = 'step', bins = 25, orientation='horizontal', log = True)
		# yhistAxis.hist(averageFoldChangePerMonomer[smallBurst], histtype = 'step', bins = 100, orientation='horizontal', log = True, color="green")

		yhistAxis.hist(averageFoldChangePerMonomer[~smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.01), np.log10(1000.), 25), range = [0.7, 100], log = True,  orientation='horizontal', color="red", linewidth=1)
		yhistAxis.hist(averageFoldChangePerMonomer[smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), range = [0.7, 100], log = True,  orientation='horizontal', color="blue", linewidth=1)
		# yhistAxis.set_ylim([0.7, 100])
		yhistAxis.set_yscale("log")


		# xhistAxis.hist(averageInitiationEventsPerMonomer[smallBurst], histtype = 'step', color = "green", bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), log = True, range = [-10., 1000.])
		# xhistAxis.hist(averageInitiationEventsPerMonomer[~smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), log = True, range = [-10., 1000.])
		# xhistAxis.set_xscale("log")

		for label in yhistAxis.xaxis.get_ticklabels()[::2]:
			label.set_visible(False)

		whitePadSparklineAxis(scatterAxis)
		whitePadSparklineAxis(yhistAxis)
		# whitePadSparklineAxis(xhistAxis)

		yhistAxis.set_yticks([1., 2.])
		yhistAxis.set_yticklabels([])

		# Label 1 and 2 with arrows
		yhistAxis.annotate("1", xy = (1e4, 1), xytext = (1e5, 1), fontsize = FONT_SIZE, arrowprops = dict(facecolor = "black", edgecolor = "none", width = 0.5, headwidth = 4),  verticalalignment = "center")
		yhistAxis.annotate("2", xy = (1e4, 2), xytext = (1e5, 2), fontsize = FONT_SIZE, arrowprops = dict(facecolor = "black", edgecolor = "none", width = 0.5, headwidth = 4),  verticalalignment = "center")

		for tick in scatterAxis.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in scatterAxis.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

		for tick in yhistAxis.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in yhistAxis.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

		scatterAxis.tick_params(
			axis='both',          # which axis
			which='both',      # both major and minor ticks are affected
			right=False,       # ticks along the bottom edge are off
			left=True,         # ticks along the top edge are off
			top = False,
			bottom = True,
			)

		yhistAxis.tick_params(
			axis='both',          # which axis
			which='both',      # both major and minor ticks are affected
			right=False,       # ticks along the bottom edge are off
			left=True,         # ticks along the top edge are off
			top = False,
			bottom = True,
			)

		# xhistAxis.tick_params(
		# 	axis='both',          # which axis
		# 	which='both',      # both major and minor ticks are affected
		# 	right=False,      # ticks along the bottom edge are off
		# 	left=True,         # ticks along the top edge are off
		# 	top=False,
		# 	bottom=True,
		# 	)

		plt.subplots_adjust(wspace=0.3, hspace=0.3, bottom = 0.2)

		# for label in yhistAxis_2.xaxis.get_ticklabels()[::2]:
		# 	label.set_visible(False)

		# axesList[0].semilogy(uniqueBurstSizes, probExistByBurstSize)
		# axesList[1].semilogy(uniqueBurstSizes, probDoubleByBurstSize)

		# # axesList[0].set_ylabel("Probability exists")
		# # axesList[1].set_ylabel("Probability doubles")
		# # axesList[1].set_xlabel("Number of transcription events per generation")

		# axesList[2].semilogy(uniqueBurstSizes, probExistByBurstSize)
		# axesList[2].set_xlim([0., 10.])
		# #axesList[2].set_ylim([0.96, 1.0])
		# axesList[3].semilogy(uniqueBurstSizes, probDoubleByBurstSize)
		# axesList[3].set_xlim([0., 10.])
		# #axesList[3].set_ylim([0.96, 1.0])

		# axesList[0].set_ylabel("Probability\nexists")
		# axesList[1].set_ylabel("Probability\ndoubles")
		# axesList[2].set_ylabel("Probability\nexists")
		# axesList[3].set_ylabel("Probability\ndoubles")
		# axesList[3].set_xlabel("Number of transcription events per generation")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		#plt.close("all")

		scatterAxis.set_xlabel("")
		scatterAxis.set_ylabel("")
		scatterAxis.set_xticklabels([])
		scatterAxis.set_yticklabels([])

		yhistAxis.set_xlabel("")
		yhistAxis.set_ylabel("")
		yhistAxis.set_xticklabels([])

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
