"""
Plot protein monomer counts

@author: John Mason, Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


# TODO: account for complexation

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get the names of proteins from the KB

		sim_data = cPickle.load(open(simDataFile, "rb"))

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

		avgCounts = view_translation.counts()

		relativeCounts = avgCounts / avgCounts.sum()

		expectedCountsArbitrary = normalize(
			sim_data.process.transcription.rnaExpression[sim_data.condition][sim_data.relation.rnaIndexToMonomerMapping] *
			sim_data.process.translation.translationEfficienciesByMonomer /
			(np.log(2) / sim_data.doubling_time.asNumber(units.s) + sim_data.process.translation.monomerData["degRate"].asNumber(1 / units.s))
			)

		expectedCountsRelative = expectedCountsArbitrary / expectedCountsArbitrary.sum()

		plt.figure(figsize = (8.5, 11))

		maxLine = 1.1 * max(np.log10(expectedCountsRelative.max() + 1), np.log10(relativeCounts.max() + 1))
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1), 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("log10(Expected protein distribution (from fitting))")
		plt.ylabel("log10(Actual protein distribution (average over life cycle))")
		plt.title("PCC (of log values): %0.2f" % pearsonr(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1))[0])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
