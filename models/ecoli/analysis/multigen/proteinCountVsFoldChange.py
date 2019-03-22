from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

FROM_CACHE = False

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all ids reqiured
		sim_data = cPickle.load(open(simDataFile, "rb"))
		ids_complexation = sim_data.process.complexation.moleculeNames # Complexe of proteins, and protein monomers
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes # Only complexes
		ids_equilibrium = sim_data.process.equilibrium.moleculeNames # Complexes of proteins + small molecules, small molecules, protein monomers
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes # Only complexes
		ids_translation = sim_data.process.translation.monomerData["id"].tolist() # Only protein monomers

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
		monomerAvgCountMultigen = np.zeros((n_sims, n_monomers), dtype=np.float)

		if not FROM_CACHE:

			for gen_idx, simDir in enumerate(allDir):
				simOutDir = os.path.join(simDir, "simOut")

				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				## READ DATA ##
				# Read in bulk ids and counts
				bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

				if first_build:
					moleculeIds = bulkMolecules.readAttribute("objectNames")

					moleculeDict = {mol: i for i, mol in enumerate(moleculeIds)}
					complexationIdx = np.array([moleculeDict[x] for x in ids_complexation]) # Complexe of proteins, and protein monomers
					complexation_complexesIdx = np.array([moleculeDict[x] for x in ids_complexation_complexes]) # Only complexes
					equilibriumIdx = np.array([moleculeDict[x] for x in ids_equilibrium]) # Complexes of proteins + small molecules, small molecules, protein monomers
					equilibrium_complexesIdx = np.array([moleculeDict[x] for x in ids_equilibrium_complexes]) # Only complexes
					translationIdx = np.array([moleculeDict[x] for x in ids_translation]) # Only protein monomers

					ribosomeIdx = np.array([moleculeDict[x] for x in ribosome_subunit_ids])
					rnapIdx = np.array([moleculeDict[x] for x in rnap_subunit_ids])

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

				bulkCounts[:, ribosomeIdx] += ribosomeSubunitCounts.astype(np.int)
				bulkCounts[:, rnapIdx] += rnapSubunitCounts.astype(np.int)

				# Get protein monomer counts for calculations now that all complexes are dissociated
				proteinMonomerCounts = bulkCounts[:, translationIdx]

				## CALCULATIONS ##
				monomerAverageCount = proteinMonomerCounts.mean(axis=0)

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

				# Log data
				monomerExistMultigen[gen_idx,:] = monomerExist
				ratioFinalToInitialCountMultigen[gen_idx,:] = ratioFinalToInitialCount
				initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer
				monomerAvgCountMultigen[gen_idx, :] = monomerAverageCount

			cPickle.dump(monomerExistMultigen, open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "wb"))
			cPickle.dump(ratioFinalToInitialCountMultigen, open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "wb"))
			cPickle.dump(initiationEventsPerMonomerMultigen, open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "wb"))
			cPickle.dump(monomerAvgCountMultigen, open(os.path.join(plotOutDir,"monomerAvgCountMultigen.pickle"), "wb"))

		monomerExistMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerExistMultigen.pickle"), "rb"))
		ratioFinalToInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "rb"))
		initiationEventsPerMonomerMultigen = cPickle.load(open(os.path.join(plotOutDir,"initiationEventsPerMonomerMultigen.pickle"), "rb"))
		monomerAvgCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"monomerAvgCountMultigen.pickle"), "rb"))


		existFractionPerMonomer = monomerExistMultigen.mean(axis=0)
		averageFoldChangePerMonomer = ratioFinalToInitialCountMultigen.mean(axis=0)
		averageInitiationEventsPerMonomer = initiationEventsPerMonomerMultigen.mean(axis=0)
		averageCountPerMonomer = monomerAvgCountMultigen.mean(axis=0)

		translationEff = sim_data.process.translation.translationEfficienciesByMonomer

		# import ipdb; ipdb.set_trace()
		# uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		# probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		# probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		# for idx, burstSize in enumerate(uniqueBurstSizes):
		# 	mask = initiationEventsPerMonomerMultigen == burstSize
		# 	mask_sum = mask.sum()
		# 	probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
		# 	probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())


		# fig, axesList = plt.subplots(4,1)

		scatterAxis = plt.subplot2grid((4,4), (1, 0), colspan=3, rowspan=3)#, sharex = xhistAxis, sharey = yhistAxis)
		xhistAxis = plt.subplot2grid((4,4), (0,0), colspan=3, sharex = scatterAxis)
		yhistAxis = plt.subplot2grid((4,4), (1,3), rowspan=3, sharey = scatterAxis)

		xhistAxis.xaxis.set_visible(False)
		yhistAxis.yaxis.set_visible(False)

		scatterAxis.loglog(averageCountPerMonomer, averageFoldChangePerMonomer, marker = '.', alpha = 0.9, lw = 0.)#, s = 5)
		scatterAxis.set_ylabel("Average fold change per monomer per generation")
		scatterAxis.set_xlabel("Average monomer count per generation per monomer")

		scatterAxis.set_ylim([10e-1, 10e1])
		scatterAxis.set_xlim([10e-3, 10e5])

		yhistAxis.hist(averageFoldChangePerMonomer, bins = np.logspace(np.log10(10e-1), np.log10(10e1), 125), orientation='horizontal', log = True, range = [10e-3, 10e6])
		yhistAxis.set_yscale("log")

		xhistAxis.hist(averageCountPerMonomer, bins = np.logspace(np.log10(10e-3), np.log10(10e5), 125), log = True, range = [1e-3, 1e5])
		xhistAxis.set_xscale("log")



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

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
