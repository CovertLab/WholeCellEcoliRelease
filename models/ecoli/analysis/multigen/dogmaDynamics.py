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

FROM_CACHE = False

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

	ids_transcription = sim_data.process.transcription.rnaData["id"].tolist()

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

	ratioFinalToInitialCountMultigen = np.zeros((n_sims, n_monomers), dtype = np.float)
	initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)

	protein_index_of_interest = np.zeros(0, dtype = np.int)

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

			ratioFinalToInitialCount = (proteinMonomerCounts[-1,:] + 1) / (proteinMonomerCounts[0,:].astype(np.float) + 1)

			index_of_interest = np.where(ratioFinalToInitialCount > 100.)
			if len(index_of_interest):
				protein_index_of_interest = np.hstack((protein_index_of_interest, index_of_interest[0]))
			
		cPickle.dump(protein_index_of_interest, open(os.path.join(plotOutDir,"protein_index_of_interest.pickle"), "wb"))


	protein_index_of_interest = cPickle.load(open(os.path.join(plotOutDir,"protein_index_of_interest.pickle"), "rb"))

	complexationIdx = cPickle.load(open(os.path.join(plotOutDir,"complexationIdx.pickle"), "rb"))
	complexation_complexesIdx = cPickle.load(open(os.path.join(plotOutDir,"complexation_complexesIdx.pickle"), "rb"))
	equilibriumIdx = cPickle.load(open(os.path.join(plotOutDir,"equilibriumIdx.pickle"), "rb"))
	equilibrium_complexesIdx = cPickle.load(open(os.path.join(plotOutDir,"equilibrium_complexesIdx.pickle"), "rb"))
	translationIdx = cPickle.load(open(os.path.join(plotOutDir,"translationIdx.pickle"), "rb"))
	transcriptionIdx = cPickle.load(open(os.path.join(plotOutDir,"transcriptionIdx.pickle"), "rb"))
	ribosomeIdx = cPickle.load(open(os.path.join(plotOutDir,"ribosomeIdx.pickle"), "rb"))
	rnapIdx = cPickle.load(open(os.path.join(plotOutDir,"rnapIdx.pickle"), "rb"))

	fig, axesList = plt.subplots(2,1, sharex = True)

	for gen_idx, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		## READ DATA ##
		# Read in bulk ids and counts
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

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
		rnaMonomerCounts = bulkCounts[:, transcriptionIdx]

		if len(protein_index_of_interest):
			axesList[0].plot(time / 60., proteinMonomerCounts[:, protein_index_of_interest[0]])

			rna_index_of_interest = sim_data.relation.monomerIndexToRnaMapping[protein_index_of_interest[0]]

			axesList[1].plot(time / 60., rnaMonomerCounts[:, rna_index_of_interest])


	axesList[1].set_xlabel("Time (min)")
	axesList[0].set_ylabel("Protein monomer count")
	axesList[1].set_ylabel("mRNA count")



		# 	# Log data

		# cPickle.dump(proteinMonomerCountsMultigen, open(os.path.join(plotOutDir,"proteinMonomerCountsMultigen.pickle"), "wb"))
		# cPickle.dump(rnaMonomerCountsMultigen, open(os.path.join(plotOutDir,"rnaMonomerCountsMultigen.pickle"), "wb"))


	# import ipdb; ipdb.set_trace()
	# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
	# ratioByBurstSize = np.zeros(uniqueBurstSizes.size)

	# burstSizeToPlot = np.zeros(0)
	# ratioToPlot = np.zeros(0)
	# for idx, burstSize in enumerate(uniqueBurstSizes):
	# 	mask = initiationEventsPerMonomerMultigen == burstSize
	# 	burstSizeToPlot = np.hstack((burstSizeToPlot, np.ones(mask.sum()) * burstSize))
	# 	ratioToPlot = np.hstack((ratioToPlot, ratioFinalToInitialCountMultigen[mask]))


	# real_values_mask = np.logical_not(np.logical_or(np.isnan(ratioToPlot), np.isinf(ratioToPlot)))
	# burstSizeToPlot = burstSizeToPlot[real_values_mask]
	# ratioToPlot = ratioToPlot[real_values_mask]

	# mean = ratioToPlot.mean()
	# std = ratioToPlot.std()

	# scatterAxis = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
	# xhistAxis = plt.subplot2grid((4,4), (0,0), colspan=3, sharex = scatterAxis)
	# yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

	# xhistAxis.xaxis.set_visible(False)
	# yhistAxis.yaxis.set_visible(False)

	# scatterAxis.scatter(burstSizeToPlot, ratioToPlot, marker = '.', s = 0.5, alpha = 0.75, lw = 0.05)
	# scatterAxis.set_ylabel("Protein monomer " + r"$\frac{count_f + 1}{count_i + 1}$" + "\n" + r"Clipped $[0, 10]$")
	# scatterAxis.set_xlabel("Transcription events per generation\n" + r"Clipped $[0, 1000]$")

	# scatterAxis.set_ylim([0., 10.])
	# scatterAxis.set_xlim([0., 1000.])

	# yhistAxis.hist(ratioToPlot, bins = 100, orientation='horizontal', range = [0., 10.], log = True)
	# xhistAxis.hist(burstSizeToPlot, bins = 100, log = True, range = [0., 1000.])

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
