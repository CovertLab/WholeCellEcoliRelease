#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

import cPickle

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils import units
FROM_CACHE = False

CLOSE_TO_DOUBLE = 0.1
FONT_SIZE = 9

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all ids reqiured
	sim_data = cPickle.load(open(simDataFile, "rb"))
	ids_complexation = sim_data.process.complexation.moleculeNames # Complexe of proteins, and protein monomers
	ids_complexation_complexes = [ids_complexation[i] for i in np.where((sim_data.process.complexation.stoichMatrix() == 1).sum(axis = 1))[0]] # Only complexes
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames # Complexes of proteins + small molecules, small molecules, protein monomers
	ids_equilibrium_complexes = [ids_equilibrium[i] for i in np.where((sim_data.process.equilibrium.stoichMatrix() == 1).sum(axis = 1))[0]] # Only complexes
	ids_translation = sim_data.process.translation.monomerData["id"].tolist() # Only protein monomers

	# ids_ribosome = 
	data_50s = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s50_fullComplex[0])
	data_30s = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s30_fullComplex[0])
	ribosome_subunit_ids = data_50s["subunitIds"].tolist() + data_30s["subunitIds"].tolist()
	ribosome_subunit_stoich = np.hstack((data_50s["subunitStoich"],data_30s["subunitStoich"]))

	data_rnap = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.rnapFull[0])
	rnap_subunit_ids = data_rnap["subunitIds"].tolist()
	rnap_subunit_stoich = data_rnap["subunitStoich"]

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	first_build = True

	# Pre-allocate variables. Rows = Generations, Cols = Monomers
	n_monomers = sim_data.process.translation.monomerData['id'].size
	n_sims = ap.n_generation

	monomerExistMultigen = np.zeros((n_sims, n_monomers), dtype = np.bool)
	ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = np.float)
	initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
	monomerCountInitialMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)
	cellMassInitialMultigen = np.zeros(n_sims, dtype = np.float)

	if not FROM_CACHE:

		for gen_idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

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
			bulkCounts[:, complexationIdx] += np.dot(sim_data.process.complexation.stoichMatrixMonomers(), bulkCounts[:, complexation_complexesIdx].transpose() * -1).transpose()

			# Dissociate protein-small molecule complexes
			bulkCounts[:, equilibriumIdx] += np.dot(sim_data.process.equilibrium.stoichMatrixMonomers(), bulkCounts[:, equilibrium_complexesIdx].transpose() * -1).transpose()

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
			initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rnaIndexToMonomerMapping]

			# Load cell mass
			cellMassInitial = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")[0]

			# Log data
			monomerExistMultigen[gen_idx,:] = monomerExist
			ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer
			monomerCountInitialMultigen[gen_idx,:] = proteinMonomerCounts[0,:]
			cellMassInitialMultigen[gen_idx] = cellMassInitial

		cPickle.dump(monomerExistMultigen, open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "wb"))
		cPickle.dump(ratioFinalToInitialCountMultigen, open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "wb"))
		cPickle.dump(initiationEventsPerMonomerMultigen, open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "wb"))
		cPickle.dump(monomerCountInitialMultigen, open(os.path.join(plotOutDir,"monomerCountInitialMultigen.pickle"), "wb"))
		cPickle.dump(cellMassInitialMultigen, open(os.path.join(plotOutDir,"cellMassInitialMultigen.pickle"), "wb"))


	monomerExistMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "rb"))
	ratioFinalToInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "rb"))
	initiationEventsPerMonomerMultigen = cPickle.load(open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "rb"))
	monomerCountInitialMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerCountInitialMultigen.pickle"), "rb"))
	cellMassInitialMultigen = cPickle.load(open(os.path.join(plotOutDir,"cellMassInitialMultigen.pickle"), "rb"))
	cellMassInitialMultigen = units.fg * cellMassInitialMultigen

	existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
	averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen.mean(axis=0)
	averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)


	mws = sim_data.getter.getMass(sim_data.process.translation.monomerData['id'])
	monomerInitialMasses = (mws * monomerCountInitialMultigen / sim_data.constants.nAvogadro)

	# np.tile(cellMassInitialMultigen.asNumber().reshape((1,10)), (n_monomers,1))

	# initialMassFractions = monomerInitialMasses.asNumber(units.fg).transpose() / np.tile(cellMassInitialMultigen.asNumber().reshape((1,10)), (n_monomers,1))
	# averageInitialMassFractions = initialMassFractions.mean(axis = 1)
	avgMonomerInitialMass = monomerInitialMasses.asNumber(units.fg).mean(axis=0)
	avgMonomerInitialMassFraction = avgMonomerInitialMass / avgMonomerInitialMass.sum()

	# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
	# probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
	# probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

	# for idx, burstSize in enumerate(uniqueBurstSizes):
	# 	mask = initiationEventsPerMonomerMultigen == burstSize
	# 	mask_sum = mask.sum()
	# 	probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
	# 	probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())


	# fig, axesList = plt.subplots(4,1)

	scatterAxis = plt.subplot2grid((4,5), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
	# scatterAxis.axhline(1.0, linewidth=0.5, color='black', linestyle="--", xmin = 0.5, xmax = 1.)
	# scatterAxis.axhline(2.0, linewidth=0.5, color='black', linestyle="--", xmin = 0.5, xmax = 1.)
	xhistAxis = plt.subplot2grid((4,5), (0,0), colspan=3, sharex = scatterAxis)
	yhistAxis = plt.subplot2grid((4,5), (1,3), rowspan=3, sharey = scatterAxis)
	yhistAxis.axhline(1.0, linewidth=0.5, color='black', linestyle="--")
	yhistAxis.axhline(2.0, linewidth=0.5, color='black', linestyle="--")
	yhistAxis_2 = plt.subplot2grid((4,5), (1,4), rowspan=3, sharey = scatterAxis)
	yhistAxis_2.axhline(1.0, linewidth=0.5, color='black', linestyle="--")
	yhistAxis_2.axhline(2.0, linewidth=0.5, color='black', linestyle="--")

	xhistAxis.xaxis.set_visible(False)
	yhistAxis.yaxis.set_visible(False)

	smallBurst = averageInitiationEventsPerMonomer <= 1.


	scatterAxis.semilogx(averageInitiationEventsPerMonomer[smallBurst], averageFoldChangePerMonomer[smallBurst], marker = '.', color = "red", alpha = 0.9, lw = 0.)#, s = 5)
	scatterAxis.semilogx(averageInitiationEventsPerMonomer[~smallBurst], averageFoldChangePerMonomer[~smallBurst], marker = '.', color = "blue", alpha = 0.9, lw = 0.)#, s = 5)
	scatterAxis.set_ylabel("Average fold change per\nmonomer per generation", fontsize = FONT_SIZE)
	scatterAxis.set_xlabel("Average number of transcription\nevents per monomer per generation", fontsize = FONT_SIZE)

	lims = yhistAxis.get_ylim()
	step = (lims[1] - lims[0]) / 125
	bins = np.arange(lims[0], lims[1] + step, step)

	mass_in_binrange = np.zeros(bins.size-1, dtype=np.float)
	for i in range(len(bins) - 1):
		in_bin_range = np.logical_and(averageFoldChangePerMonomer > bins[i], averageFoldChangePerMonomer < bins[i+1])
		mass_in_binrange[i] = avgMonomerInitialMassFraction[in_bin_range].sum()

	yhistAxis_2.barh(bottom = bins[:-1], width = mass_in_binrange, height=(lims[1] - lims[0]) / 125, color = "white")
	yhistAxis_2.set_xlim([0., 1.])
	yhistAxis_2.yaxis.set_label_position("right")
	yhistAxis_2.set_xlabel("Fraction of\nproteome mass", fontsize = FONT_SIZE)
	scatterAxis.set_xlim([-10., 1000.])

	yhistAxis.hist(averageFoldChangePerMonomer[~smallBurst], histtype = 'step', bins = 25, orientation='horizontal', log = True)
	yhistAxis.hist(averageFoldChangePerMonomer[smallBurst], histtype = 'step', bins = 100, orientation='horizontal', log = True, color="red")
	yhistAxis_2.yaxis.tick_right()

	xhistAxis.hist(averageInitiationEventsPerMonomer[smallBurst], histtype = 'step', color = "red", bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), log = True, range = [-10., 1000.])
	xhistAxis.hist(averageInitiationEventsPerMonomer[~smallBurst], histtype = 'step', bins = np.logspace(np.log10(0.01), np.log10(1000.), 125), log = True, range = [-10., 1000.])
	xhistAxis.set_xscale("log")

	for label in yhistAxis.xaxis.get_ticklabels()[::2]:
		label.set_visible(False)

	for label in yhistAxis_2.xaxis.get_ticklabels()[::2]:
		label.set_visible(False)

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


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "fixme")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
