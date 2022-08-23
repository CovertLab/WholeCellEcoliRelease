"""
Plots average RNA counts per cell vs average protein counts per cell.
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
import matplotlib.lines as mlines
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

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

PLOT_ZEROS_ON_LINE = 2.5e-6


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		HIGHLIGHT_GENES = False

		# Check if cache from figure5B_E_F_G.py exist
		if os.path.exists(os.path.join(plotOutDir, "figure5B.pickle")):
			figure5B_data = cPickle.load(open(os.path.join(plotOutDir, "figure5B.pickle"), "rb"))
			colors = figure5B_data["colors"]
			mrnaIds = figure5B_data["id"].tolist()
		else:
			print("Requires figure5B.pickle from figure5B_E_F_G.py")
			return

		# Get all cells
		allDir = self.ap.get_cells()

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rna_data["id"][sim_data.relation.cistron_to_monomer_mapping] # orders rna IDs to match monomer IDs

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
		monomersInvolvedInManyComplexes = []
		monomersInvolvedInComplexes = []
		for complexId in ids_complexation_complexes:
			subunitIds = sim_data.process.complexation.get_monomers(complexId)["subunitIds"]
			for subunitId in subunitIds:
				if subunitId in monomersInvolvedInComplexes:
					monomersInvolvedInManyComplexes.append(subunitId)
				monomersInvolvedInComplexes.append(subunitId)
		monomersInvolvedInManyComplexes_id = list(set(monomersInvolvedInManyComplexes))
		monomersInvolvedInManyComplexes_dict = {}
		for x in monomersInvolvedInManyComplexes_id:
			monomersInvolvedInManyComplexes_dict[x] = {}

		# Get average (over timesteps) counts for All genseration (ie. All cells)
		avgRnaCounts_forAllCells = np.zeros(rnaIds.shape[0], np.float64)
		avgProteinCounts_forAllCells = np.zeros(rnaIds.shape[0], np.float64)
		for i, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			# Account for bulk molecules
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeCounts = bulkMolecules.readColumn("counts")
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], int)
			proteinCountsBulk = bulkMoleculeCounts[:, proteinIndexes]
			rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], int)
			avgRnaCounts = bulkMoleculeCounts[:, rnaIndexes].mean(axis = 0)
			bulkMolecules.close()
			if i == 0:
				# Skip first few time steps for 1st generation (becaused complexes have not yet formed during these steps)
				bulkContainer.countsIs(np.mean(proteinCountsBulk[5:, :], axis = 0))
			else:
				bulkContainer.countsIs(proteinCountsBulk.mean(axis = 0))

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

			# Average counts of monomers
			avgMonomerCounts = view_translation.counts()

			# Get counts of "functional units" (ie. complexed forms)
			avgProteinCounts = avgMonomerCounts[:]
			avgComplexCounts = view_complexation_complexes.counts()

			for j, complexId in enumerate(ids_complexation_complexes):
				# Map all subsunits to the average counts of the complex (ignores counts of monomers)
				# Some subunits are involved in multiple complexes - these cases are kept track
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

					if subunitId not in monomersInvolvedInManyComplexes_id:
						avgProteinCounts[ids_translation.index(subunitId)] = avgComplexCounts[j]
					else:
						if complexId not in monomersInvolvedInManyComplexes_dict[subunitId]:
							monomersInvolvedInManyComplexes_dict[subunitId][complexId] = 0.
						monomersInvolvedInManyComplexes_dict[subunitId][complexId] += avgComplexCounts[j]

			# Store
			avgRnaCounts_forAllCells += avgRnaCounts
			avgProteinCounts_forAllCells += avgProteinCounts


		# Per cell
		avgRnaCounts_perCell = avgRnaCounts_forAllCells / float(len(allDir))
		avgProteinCounts_perCell = avgProteinCounts_forAllCells / float(len(allDir))

		# Plot
		fig, ax = plt.subplots(1, 1, figsize = (10, 10))

		for monomer in monomersInvolvedInManyComplexes_id:
			index = ids_translation.index(monomer)
			color_index = mrnaIds.index(rnaIds[index])
			color = colors[color_index]

			for complexId in monomersInvolvedInManyComplexes_dict[monomer]:
				avgComplexCount = monomersInvolvedInManyComplexes_dict[monomer][complexId] / float(len(allDir))

				if avgComplexCount == 0:
					ax.loglog(avgRnaCounts_perCell[index], 2.5e-6, alpha = 0.5, marker = ".", lw = 0., color = color)

				else:
					if avgRnaCounts_perCell[index] == 0:
						ax.loglog(PLOT_ZEROS_ON_LINE, avgComplexCount, alpha = 0.5, marker = ".", lw = 0., color = color)
					else:
						ax.loglog(avgRnaCounts_perCell[index], avgComplexCount, alpha = 0.5, marker = ".", lw = 0., color = color)

		# plot monomers that are not involved in complexes or involved in only 1 complex
		monomersInvolvedInManyComplexes_index = [ids_translation.index(x) for x in monomersInvolvedInManyComplexes_id]
		A = [x for x in range(len(ids_translation)) if x not in monomersInvolvedInManyComplexes_index]
		for i in A:
			color = colors[mrnaIds.index(rnaIds[i])]
			ax.loglog(avgRnaCounts_perCell[i], avgProteinCounts_perCell[i], alpha = 0.5, marker = ".", lw = 0., color = color)
		# ax.loglog(avgRnaCounts_perCell[A], avgProteinCounts_perCell[A], alpha = 0.5, marker = ".", lw = 0., color = plot_colors)

		# Plot genes with zero transcripts an arbitrary line
		noTranscripts_indices = [x for x in np.where(avgRnaCounts_perCell == 0)[0] if x not in monomersInvolvedInManyComplexes_index]
		for i in noTranscripts_indices:
			color = colors[mrnaIds.index(rnaIds[i])]
			ax.loglog(PLOT_ZEROS_ON_LINE, avgProteinCounts_perCell[i], alpha = 0.5, marker = ".", lw = 0., color = color)

		# Highlight
		if HIGHLIGHT_GENES:
			rnaIds = rnaIds.tolist()
			highlights_rnaId = ["EG12437_RNA[c]", "EG12058_RNA[c]"] # menE, ccmB
			colors = ["g", "r"]
			for i, rna in enumerate(highlights_rnaId):
				if avgRnaCounts_perCell[rnaIds.index(rna)] == 0:
					ax.loglog(PLOT_ZEROS_ON_LINE, avgProteinCounts_perCell[rnaIds.index(rna)], marker = '.', lw = 0., color = colors[i], ms = 15)
				else:
					ax.loglog(avgRnaCounts_perCell[rnaIds.index(rna)], avgProteinCounts_perCell[rnaIds.index(rna)], marker = '.', lw = 0., color = colors[i], ms = 15)

			green_dot = mlines.Line2D([], [], color = "green", linewidth = 0., marker = ".", markersize = 15, label = "menE")
			red_dot = mlines.Line2D([], [], color = "red", linewidth = 0., marker = ".", markersize = 15, label = "ccmB")
			plt.legend(handles = [green_dot, red_dot], loc = "lower right")

		# ax.hlines(1, ax.get_xlim()[0], ax.get_xlim()[1], linestyle = "--")
		ax.hlines(9786.77, ax.get_xlim()[0], ax.get_xlim()[1], linestyle = "--")

		ax.set_title("Each (translatable) gene's functional unit is represented as a point\n(ie. x points per gene where x == number of complexes the monomer is involved in)\n(avg across %s generations)" % len(allDir))
		ax.set_xlabel("<RNA> per cell")
		ax.set_ylabel("<Functional units (protein)> per cell")
		ax.tick_params(which = "both", direction = "out")

		plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.9, right = 0.95)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
