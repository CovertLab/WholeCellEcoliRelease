"""
Compare protein counts to Wisniewski 2014 data set

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/3/2015
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
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

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
		view_validation = bulkContainer.countsView(validation_data.protein.wisniewski2014Data["monomerId"].tolist())
		view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())

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

		wisniewskiCounts = validation_data.protein.wisniewski2014Data["avgCounts"]

		fig, ax = plt.subplots(2, sharey=True, figsize = (8.5, 11))

		# Wisniewski Counts
		ax[0].scatter(np.log10(wisniewskiCounts + 1), np.log10(view_validation.counts() + 1), c='w', edgecolor = 'k', alpha=.7)
		ax[0].set_xlabel("log10(Wisniewski 2014 Counts)")
		ax[0].set_title("Pearson r: %0.2f" % pearsonr(np.log10(view_validation.counts() + 1), np.log10(wisniewskiCounts + 1))[0])

		# Schmidt Counts
		schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]
		ax[1].scatter(
			np.log10(schmidtCounts + 1),
			np.log10(view_validation_schmidt.counts() + 1),
			c='w', edgecolor = 'k', alpha=.7)
		ax[1].set_xlabel("log10(Schmidt 2015 Counts)")
		ax[1].set_title("Pearson r: %0.2f" % pearsonr(np.log10(view_validation_schmidt.counts() + 1), np.log10(schmidtCounts + 1))[0])

		plt.ylabel("log10(Simulation Average Counts)")
		# NOTE: This Pearson correlation goes up (at the time of writing) about 0.05 if you only
		# include proteins that you have translational efficiencies for
		plt.xlim(xmin=0)
		plt.ylim(ymin=0)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
