from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt
from six.moves import cPickle, range

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

PLOT_ZEROS_ON_LINE = 2.5e-6

complexToMonomer = {
	"CPLX0-7620[c]": "PD00260[c]", # CPLX0-7620's monomer is EG10359-MONOMER, which is ID'ed as PD00260 (proteins.tsv)
	"CPLX0-8801[c]": "G6420-MONOMER[c]",
	"CPLX0-7677[c]": "EG11969-MONOMER[c]",
	"CPLX0-7702[c]": "G6263-MONOMER[c]",
	"CPLX0-7701[c]": "G6420-MONOMER[c]",
	}

monomerToTranslationMonomer = {
	"MONOMER0-1781[c]": "EG11519-MONOMER[c]", # MONOMER0-1781 is a complex, EG11519 is its monomer
	"EG10359-MONOMER[c]": "PD00260[c]", # EG10359-MONOMER is not the ID of fur monomer, it's PD00260 (proteins.tsv)
	}


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		# Check if cache from rnaVsProteinPerCell.py exists
		if os.path.exists(os.path.join(plotOutDir, "rnaVsProteinPerCell_alltimesteps.cPickle")):
			rnaVsProteinPerCell = cPickle.load(open(os.path.join(plotOutDir, "rnaVsProteinPerCell_alltimesteps.cPickle"), "rb"))
			avgProteinCounts_forAllCells = rnaVsProteinPerCell["protein"]
			monomersInManyComplexes_avgProteins_dict = rnaVsProteinPerCell["monomersInManyComplexes"]

			avgProteinCounts_perCell = avgProteinCounts_forAllCells / float(32)
			# Copy structure of monomersInManyComplexes_dict
			monomersInManyComplexes_dict = monomersInManyComplexes_avgProteins_dict.copy()

			for key in monomersInManyComplexes_dict.keys():
				monomersInManyComplexes_dict[key] = {}
		else:
			print("Requires rnaVsProteinPerCell.cPickle from rnaVsProteinPerCell.py")
			return

		# Check if cache from figure5B_E_F_G.py exist
		if os.path.exists(os.path.join(plotOutDir, "figure5B.pickle")):
			figure5B_data = cPickle.load(open(os.path.join(plotOutDir, "figure5B.pickle"), "rb"))
			colors = figure5B_data["colors"]
			mrnaIds = figure5B_data["id"].tolist()
		else:
			print("Requires figure5B.pickle from figure5B_E_F_G.py")
			return

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rna_data["id"][sim_data.relation.RNA_to_monomer_mapping] # orders rna IDs to match monomer IDs

		# Make views for monomers
		ids_complexation = sim_data.process.complexation.molecule_names
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes
		ids_equilibrium = sim_data.process.equilibrium.molecule_names
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes
		ids_translation = sim_data.process.translation.monomer_data["id"].tolist()
		ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))
		bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
		view_complexation = bulkContainer.countsView(ids_complexation)
		view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
		view_equilibrium = bulkContainer.countsView(ids_equilibrium)
		view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
		view_translation = bulkContainer.countsView(ids_translation)

		# Identify monomers that are subunits for multiple complexes
		monomersInManyComplexes = []
		monomersInComplexes = []
		for complexId in ids_complexation_complexes:
			subunitIds = sim_data.process.complexation.get_monomers(complexId)["subunitIds"]
			for subunitId in subunitIds:
				if subunitId in monomersInComplexes:
					monomersInManyComplexes.append(subunitId)
				monomersInComplexes.append(subunitId)
		monomersInManyComplexes_id = list(set(monomersInManyComplexes))

		# Initialize minimum number (counts) of each functional unit
		minProteinCounts = np.ones(rnaIds.shape[0], np.float64) * np.inf

		for i, simDir in enumerate(allDir):
			print(i)
			simOutDir = os.path.join(simDir, "simOut")

			# Account for bulk molecules
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], np.int)
			proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]
			bulkMolecules.close()

			if i == 0:
				# Skip first few time steps for 1st generation (becaused complexes have not yet formed during these steps)
				bulkContainer.countsIs(np.min(proteinCountsBulk[5:, :], axis = 0))
			else:
				bulkContainer.countsIs(proteinCountsBulk.min(axis = 0))

			# Unique molecules
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
			rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_RNAP')
			nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
			uniqueMoleculeCounts.close()
			# Account for unique molecules
			bulkContainer.countsInc(nActiveRibosome.mean(), [sim_data.molecule_ids.s30_full_complex, sim_data.molecule_ids.s50_full_complex])
			bulkContainer.countsInc(nActiveRnaPoly.mean(), [sim_data.molecule_ids.full_RNAP])
			# Account for small-molecule bound complexes
			view_equilibrium.countsInc(np.dot(sim_data.process.equilibrium.stoich_matrix_monomers(), view_equilibrium_complexes.counts() * -1))

			# Minimum monomer counts
			minMonomerCounts = view_translation.counts()

			# Get counts of "functional units" (ie. complexed forms)
			minProteinCounts_thisGen = minMonomerCounts[:]
			minComplexCounts = view_complexation_complexes.counts()

			for j, complexId in enumerate(ids_complexation_complexes):
				# Map all subsunits to the minimum counts of the complex (ignores counts of monomers)
				# Some subunits are involved in multiple complexes - these cases are kept track separately
				# by monomersInManyComplexes_dict
				subunitIds = sim_data.process.complexation.get_monomers(complexId)["subunitIds"]

				for subunitId in subunitIds:
					if subunitId not in ids_translation:
						if subunitId in monomerToTranslationMonomer:
							# couple monomers have different ID in ids_translation
							subunitId = monomerToTranslationMonomer[subunitId]
						elif "CPLX" in subunitId:
							# few transcription factors are complexed with ions
							subunitId = complexToMonomer[subunitId]
						elif "RNA" in subunitId:
							continue

					if subunitId not in monomersInManyComplexes_id:
						minProteinCounts_thisGen[ids_translation.index(subunitId)] = minComplexCounts[j]
					else:
						minProteinCounts_thisGen[ids_translation.index(subunitId)] = np.inf

						if complexId not in monomersInManyComplexes_dict[subunitId]:
							monomersInManyComplexes_dict[subunitId][complexId] = np.inf

						prev_entry = monomersInManyComplexes_dict[subunitId][complexId]
						monomersInManyComplexes_dict[subunitId][complexId] = min(prev_entry, minComplexCounts[j])

			# Store
			minProteinCounts = np.minimum(minProteinCounts, minProteinCounts_thisGen)

		# Plot
		plt.figure(figsize = (11, 8.5))
		ax0 = plt.subplot2grid((3, 2), (0, 0))
		ax2 = plt.subplot2grid((3, 2), (1, 0))
		ax3 = plt.subplot2grid((3, 2), (0, 1))
		ax4 = plt.subplot2grid((3, 2), (1, 1))
		ax5 = plt.subplot2grid((3, 2), (2, 1))
		axesList = [ax0, ax2, ax3, ax4, ax5]

		# Plot
		zero_functionalUnits = []
		for monomer in monomersInManyComplexes_id:
			index = ids_translation.index(monomer)
			color = colors[mrnaIds.index(rnaIds[index])]

			for complexId in monomersInManyComplexes_dict[monomer]:
				minComplexCount = monomersInManyComplexes_dict[monomer][complexId]

				if minComplexCount == 0:
					zero_functionalUnits.append(avgProteinCounts_perCell[index])
					minComplexCount = 0.1 # plot 0 on arbitrary line

				avgProtein = monomersInManyComplexes_avgProteins_dict[monomer][complexId]
				if avgProtein == 0.:
					avgProtein = 0.1 # plot 0 on arbitrary line
				ax0.loglog(avgProtein, minComplexCount, alpha = 0.5, color = color, lw = 0., marker = ".")

				# plot each color (ie. subgeneration category) separately
				if color == "y":
					ax3.loglog(avgProtein, minComplexCount, color = color, markeredgecolor = "k", lw = 0., marker = ".")
				elif color == "b":
					ax4.loglog(avgProtein, minComplexCount, color = color, markeredgecolor = "k", lw = 0., marker = ".")
				else:
					ax5.loglog(avgProtein, minComplexCount, color = color, markeredgecolor = "k", lw = 0., marker = ".")

		# plot monomers that are not involved in complexes or involved in only 1 complex
		monomersInManyComplexes_index = [ids_translation.index(x) for x in monomersInManyComplexes_id]
		A = [x for x in range(len(ids_translation)) if x not in monomersInManyComplexes_index]
		for i in A:
			color = colors[mrnaIds.index(rnaIds[i])]
			ax0.loglog(avgProteinCounts_perCell[i], minProteinCounts[i], alpha = 0.5, color = color,lw = 0., marker = ".")

			if color == "y":
				ax3.loglog(avgProteinCounts_perCell[i], minProteinCounts[i], color = color, lw = 0., marker = ".")
			elif color == "b":
				ax4.loglog(avgProteinCounts_perCell[i], minProteinCounts[i], color = color, lw = 0., marker = ".")
			else:
				ax5.loglog(avgProteinCounts_perCell[i], minProteinCounts[i], color = color, lw = 0., marker = ".")

		zero_functionalUnits += avgProteinCounts_perCell[A][np.where(minProteinCounts[A] == 0)[0]].tolist()

		noProteinIndex = [x for x in np.where(minProteinCounts == 0)[0] if x not in monomersInManyComplexes_index]
		for i in noProteinIndex:
			color = colors[mrnaIds.index(rnaIds[i])]
			ax0.loglog(avgProteinCounts_perCell[i], 0.1 , alpha = 0.5, color = color,lw = 0., marker = ".")

			if color == "y":
				ax3.loglog(avgProteinCounts_perCell[i], 0.1, color = color, lw = 0., marker = ".")
			elif color == "b":
				ax4.loglog(avgProteinCounts_perCell[i], 0.1, color = color, lw = 0., marker = ".")
			else:
				ax5.loglog(avgProteinCounts_perCell[i], 0.1, color = color, lw = 0., marker = ".")

		# histograms
		bins = np.logspace(-3, 6, 75)
		n2, bins2, patches2 = ax2.hist(zero_functionalUnits, bins = bins)

		ax2.set_xscale("log")
		ax2.set_ylabel("Frequency")

		ax0.set_ylim([0.025, ax0.get_ylim()[1]])
		ax0.vlines(1, ax0.get_ylim()[0], ax0.get_ylim()[1], linestyle = "--")

		for ax in [ax3, ax4, ax5]:
			ax.set_ylim(ax0.get_ylim())
			ax.set_xlim(ax0.get_xlim())
			ax.vlines(1, ax.get_ylim()[0], ax.get_ylim()[1], linestyle = "--")

		ax0.set_ylabel("min(counts of functional units)")
		ax0.set_xlabel("<counts of functional unit> per cell")
		ax2.set_xlabel("<counts of functional unit> per cell == 0")

		for ax in axesList:
			ax.tick_params(which = "both", direction = "out")

		plt.subplots_adjust(hspace = 0.3, wspace = 0.25, left = 0.1, bottom = 0.1, top = 0.95, right = 0.95)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
