from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure

import cPickle
from models.ecoli.analysis import multigenAnalysisPlot

CLOSE_TO_DOUBLE = 0.1


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		print "DISABLED"
		return

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

		# ids_ribosome =
		data_50s = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)
		data_30s = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)
		ribosome_subunit_ids = data_50s["subunitIds"].tolist() + data_30s["subunitIds"].tolist()
		ribosome_subunit_stoich = np.hstack((data_50s["subunitStoich"],data_30s["subunitStoich"]))

		data_rnap = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.rnapFull)
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
		monomerDoubleMultigen = np.zeros((n_sims, n_monomers), dtype = np.bool)
		initiationEventsPerMonomerMultigen = np.zeros((n_sims, n_monomers), dtype = np.int)

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
			ratioFinalToInitialCount = proteinMonomerCounts[-1:] / proteinMonomerCounts[0,:].astype(np.float)
			monomerDouble = ratioFinalToInitialCount > (1 - CLOSE_TO_DOUBLE)

			# Load transcription initiation event data
			rnapData = TableReader(os.path.join(simOutDir, "RnapData"))
			initiationEventsPerRna = rnapData.readColumn("rnaInitEvent").sum(axis = 0)

			# Map transcription initiation events to monomers
			initiationEventsPerMonomer = initiationEventsPerRna[sim_data.relation.rnaIndexToMonomerMapping]

			# Log data
			monomerExistMultigen[gen_idx,:] = monomerExist
			monomerDoubleMultigen[gen_idx,:] = monomerDouble
			initiationEventsPerMonomerMultigen[gen_idx,:] = initiationEventsPerMonomer

		uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen)
		probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		for idx, burstSize in enumerate(uniqueBurstSizes):
			mask = initiationEventsPerMonomerMultigen == burstSize
			mask_sum = mask.sum()
			probExistByBurstSize[idx] = monomerExistMultigen[mask].sum() / float(mask.sum())
			probDoubleByBurstSize[idx] = monomerDoubleMultigen[mask].sum() / float(mask.sum())

		# Calculate generational standard deviation (row is generation, col is burst size)
		probExistByBurstSizeGen = np.zeros((ap.n_generation, uniqueBurstSizes.size))
		probDoubleByBurstSizeGen = np.zeros((ap.n_generation, uniqueBurstSizes.size))

		for gen_idx in range(ap.n_generation):
			for idx, burstSize in enumerate(uniqueBurstSizes):
				mask = initiationEventsPerMonomerMultigen[gen_idx,:] == burstSize

				probExistByBurstSize[gen_idx, idx] = monomerExistMultigen[gen_idx,:][mask].sum() / float(mask.sum())

				probDoubleByBurstSize[gen_idx, idx] = monomerDoubleMultigen[gen_idx,:][mask].sum() / float(mask.sum())

		fig, axesList = plt.subplots(4,1)

		axesList[0].semilogy(uniqueBurstSizes, probExistByBurstSize)
		axesList[1].semilogy(uniqueBurstSizes, probDoubleByBurstSize)

		# axesList[0].set_ylabel("Probability exists")
		# axesList[1].set_ylabel("Probability doubles")
		# axesList[1].set_xlabel("Number of transcription events per generation")

		axesList[2].semilogy(uniqueBurstSizes, probExistByBurstSize)
		axesList[2].set_xlim([0., 10.])
		#axesList[2].set_ylim([0.96, 1.0])
		axesList[3].semilogy(uniqueBurstSizes, probDoubleByBurstSize)
		axesList[3].set_xlim([0., 10.])
		#axesList[3].set_ylim([0.96, 1.0])

		axesList[0].set_ylabel("Probability\nexists")
		axesList[1].set_ylabel("Probability\ndoubles")
		axesList[2].set_ylabel("Probability\nexists")
		axesList[3].set_ylabel("Probability\ndoubles")
		axesList[3].set_xlabel("Number of transcription events per generation")


		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")

		# Ignore first 5 generations

		fig, axesList = plt.subplots(4,1)

		uniqueBurstSizes = np.unique(initiationEventsPerMonomerMultigen[5:, :])
		probExistByBurstSize = np.zeros(uniqueBurstSizes.size)
		probDoubleByBurstSize = np.zeros(uniqueBurstSizes.size)

		for idx, burstSize in enumerate(uniqueBurstSizes):
			mask = initiationEventsPerMonomerMultigen[5:, :] == burstSize
			mask_sum = mask.sum()
			probExistByBurstSize[idx] = monomerExistMultigen[5:, :][mask].sum() / float(mask.sum())
			probDoubleByBurstSize[idx] = monomerDoubleMultigen[5:, :][mask].sum() / float(mask.sum())

		axesList[0].plot(uniqueBurstSizes, probExistByBurstSize)
		axesList[1].plot(uniqueBurstSizes, probDoubleByBurstSize)

		# axesList[0].set_ylabel("Probability exists")
		# axesList[1].set_ylabel("Probability doubles")
		# axesList[1].set_xlabel("Number of transcription events per generation")

		axesList[2].plot(uniqueBurstSizes, probExistByBurstSize)
		axesList[2].set_xlim([0., 10.])
		# axesList[2].set_ylim([0.96, 1.0])
		axesList[3].plot(uniqueBurstSizes, probDoubleByBurstSize)
		axesList[3].set_xlim([0., 10.])
		# axesList[3].set_ylim([0.96, 1.0])

		axesList[0].set_ylabel("Probability\nexists")
		axesList[1].set_ylabel("Probability\ndoubles")
		axesList[2].set_ylabel("Probability\nexists")
		axesList[3].set_ylabel("Probability\ndoubles")
		axesList[3].set_xlabel("Number of transcription events per generation")

		exportFigure(plt, plotOutDir, plotOutFileName + "_skip_5_gen",metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
