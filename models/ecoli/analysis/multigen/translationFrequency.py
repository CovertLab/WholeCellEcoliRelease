#!/usr/bin/env python
"""
Plots frequency of observing at least 1 protein.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/11/2017
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()
	
	sim_data = cPickle.load(open(simDataFile, "rb"))
	tcsComplexToMonomers = sim_data.process.two_component_system.complexToMonomer
	ids_complexation = sim_data.process.complexation.moleculeNames
	ids_complexation_complexes = sim_data.process.complexation.ids_complexes
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes
	ids_twoComponent = sim_data.process.two_component_system.moleculeNames.tolist()
	ids_twoComponent_complexes = sim_data.process.two_component_system.complexToMonomer.keys()
	ids_translation = sim_data.process.translation.monomerData["id"].tolist()
	ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_twoComponent + ids_translation))

	bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
	view_complexation = bulkContainer.countsView(ids_complexation)
	view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
	view_equilibrium = bulkContainer.countsView(ids_equilibrium)
	view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
	view_twoComponent = bulkContainer.countsView(ids_twoComponent)
	view_twoComponent_complexes = bulkContainer.countsView(ids_twoComponent_complexes)
	view_translation = bulkContainer.countsView(ids_translation)

	proteinPresence = []
	for simDir in allDir:
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
		bulkContainer.countsInc(nActiveRibosome.mean(), sim_data.moleculeGroups.s30_fullComplex + sim_data.moleculeGroups.s50_fullComplex)
		bulkContainer.countsInc(nActiveRnaPoly.mean(), sim_data.moleculeGroups.rnapFull)

		# Account for two-component complexes
		view_twoComponent.countsInc(
			np.dot(sim_data.process.two_component_system.stoichMatrixMonomers(), view_twoComponent_complexes.counts() * -1)
			)

		# Account for small-molecule bound complexes
		view_equilibrium.countsInc(
			np.dot(sim_data.process.equilibrium.stoichMatrixMonomers(), view_equilibrium_complexes.counts() * -1)
			)

		# Account for monomers in complexed form
		view_complexation.countsInc(
			np.dot(sim_data.process.complexation.stoichMatrixMonomers(), view_complexation_complexes.counts() * -1)
			)

		# Get boolean protein presence
		proteinCounts = view_translation.counts()
		proteinPresence.append(proteinCounts != 0)

		# Clear counts
		bulkContainer.countsIs(0)

	proteinPresence = np.array(proteinPresence)

	# Plot
	fig = plt.figure(figsize = (12, 12))
	ax = plt.subplot(1, 1, 1)
	nGens = len(allDir)
	ax.hist(np.mean(proteinPresence, axis = 0), nGens)
	ax.set_xlabel("Frequency of observing at least 1 protein copy in 1 generation", fontsize = 14)
	ax.set_ylabel("Number of proteins", fontsize = 14)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])