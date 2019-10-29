"""
Compare protein counts to Schmidt 2015 data set

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/4/2017
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

LOW_COUNT_THRESHOLD = 30

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile, "rb"))
		validation_data = cPickle.load(open(validationDataFile, "rb"))

		ids_complexation = sim_data.process.complexation.moleculeNames
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes
		ids_equilibrium = sim_data.process.equilibrium.moleculeNames
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes
		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))
		bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
		view_complexation = bulkContainer.countsView(ids_complexation)
		view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
		view_equilibrium = bulkContainer.countsView(ids_equilibrium)
		view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
		view_translation = bulkContainer.countsView(ids_translation)
		view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDir = ap.get_cells()

		View_Validation_Schmidt = []

		fig = plt.figure(figsize = (4, 4))

		for simDir in allDir:
			# print simDir

			simOutDir = os.path.join(simDir, "simOut")

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], np.int)
			proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]
			bulkMolecules.close()

			# Account for monomers
			bulkContainer.countsIs(proteinCountsBulk.mean(axis = 0))

			# Account for unique molecules
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
			nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
			uniqueMoleculeCounts.close()
			bulkContainer.countsInc(nActiveRibosome.mean(), [sim_data.moleculeIds.s30_fullComplex, sim_data.moleculeIds.s50_fullComplex])
			bulkContainer.countsInc(nActiveRnaPoly.mean(), [sim_data.moleculeIds.rnapFull])

			# Account for small-molecule bound complexes
			view_equilibrium.countsInc(
				np.dot(sim_data.process.equilibrium.stoichMatrixMonomers(), view_equilibrium_complexes.counts() * -1)
				)

			# Account for monomers in complexed form
			view_complexation.countsInc(
				np.dot(sim_data.process.complexation.stoichMatrixMonomers(), view_complexation_complexes.counts() * -1)
				)

			view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())
			View_Validation_Schmidt.append(view_validation_schmidt.counts())

		simulation_counts = (np.array(View_Validation_Schmidt)).mean(axis = 0)

		# Schmidt Counts
		schmidtLabels = validation_data.protein.schmidt2015Data["monomerId"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		# Set up mask for proteins with low counts
		low_count_mask = schmidt_counts < LOW_COUNT_THRESHOLD
		n_low_count = low_count_mask.sum()
		n_high_count = schmidt_counts.size - n_low_count
		
		# Take logs
		schmidt_counts_log = np.log10(schmidt_counts + 1)
		simulation_counts_log = np.log10(simulation_counts + 1)

		# Compute deviations
		deviation_log = np.log10(np.abs(simulation_counts - schmidt_counts))

		axis = plt.subplot(1,1,1)

		axis.plot(schmidt_counts_log, simulation_counts_log, 'o', color = "black", markersize = 6, alpha = 0.1, zorder = 1, markeredgewidth = 0.0)
		print("R^2 (all proteins) = %.3f (n = %d)" % (
			(pearsonr(simulation_counts_log, schmidt_counts_log)[0])**2,
			schmidt_counts.size
			))
		print("R^2 (low-abundance proteins) = %.3f (n = %d)" % (
			(pearsonr(simulation_counts_log[low_count_mask],
				schmidt_counts_log[low_count_mask])[0])**2,
			n_low_count
			))
		print("R^2 (high-abundance proteins) = %.3f (n = %d)" % (
			(pearsonr(simulation_counts_log[~low_count_mask],
				schmidt_counts_log[~low_count_mask])[0])**2,
			n_high_count
			))
		
		print("Average log deviation (low-abundance proteins) = %.3f" % (
			deviation_log[low_count_mask].mean()))
		print("Average log deviation (high-abundance proteins) = %.3f" % (
			deviation_log[~low_count_mask].mean()))

		maxLine = np.ceil(
			max(schmidt_counts_log.max(), simulation_counts_log.max())
			)
		plt.plot([0, maxLine], [0, maxLine], '-k')

		plt.xlim(xmin=0, xmax=maxLine)
		plt.ylim(ymin=0, ymax=maxLine)

		axis.spines["right"].set_visible(False)
		axis.spines["top"].set_visible(False)
		axis.spines["left"].set_position(("outward", 10))
		axis.spines["bottom"].set_position(("outward", 10))
		axis.tick_params(right=False)
		axis.tick_params(top=False)
		axis.tick_params(which="both", direction="out")

		axis.set_xlim([-0.07, maxLine])
		axis.set_ylim([-0.07, maxLine])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		axis.set_yticklabels([])
		axis.set_ylabel("")
		axis.set_xticklabels([])
		axis.set_xlabel("")

		exportFigure(plt, plotOutDir, plotOutFileName + "_clean", metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
