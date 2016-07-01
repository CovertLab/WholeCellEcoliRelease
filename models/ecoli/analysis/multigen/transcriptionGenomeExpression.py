#!/usr/bin/env python
"""
Plots degree of mRNAs transcribed (out of all genes to be transcribed) for all generations.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/30/2016
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

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	# Get IDs of mRNAs
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	basalExpression = sim_data.process.transcription.rnaExpression["basal"]
	mRnaIds = np.where(isMRna)[0]

	mRnaBasalExpression = np.array([basalExpression[x] for x in mRnaIds])
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# Sort in order of decreasing basal expression
	descendingOrderIndexing = np.argsort(mRnaBasalExpression)[::-1]
	mRnaBasalExpressionSorted = mRnaBasalExpression[descendingOrderIndexing]
	mRnaNamesSorted = mRnaNames[descendingOrderIndexing]

	# Get proteins of interest
	dcuRId = np.where([x == "G7826_RNA[c]" for x in mRnaNamesSorted])[0]
	baeRId = np.where([x == "EG11618_RNA[c]" for x in mRnaNamesSorted])[0]
	narLId = np.where([x == "EG10643_RNA[c]" for x in mRnaNamesSorted])[0]

	proteinsOfInterest = np.array([dcuRId, baeRId, narLId])
	proteinsOfInterestNames = np.array(["dcuR", "baeR", "narL"])
	
	# Get expression level of mRNAs
	mRnaExpression = []

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNamesSorted])
		moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]
		bulkMolecules.close()

		moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
		mRnaExpression.append(moleculeCountsSumOverTime)

	mRnaExpression = np.array(mRnaExpression)

	# Plot
	numGens = allDir.shape[0]
	numMRnas = mRnaNamesSorted.shape[0]
	rows = numGens
	cols = 1
	fig = plt.figure(figsize = (15, 18))

	for gen in np.arange(numGens):
		ax = plt.subplot(rows, cols, gen + 1)
		ax.scatter(np.arange(numMRnas), mRnaExpression[gen], facecolor = "none", edgecolor = "b")

		heighOffset = np.max(mRnaExpression[gen]) * 0.01

		for proteinIndex, protein in enumerate(proteinsOfInterest):
			ax.scatter(protein, mRnaExpression[gen][protein], facecolor = "orange", edgecolor = "none")
			ax.annotate(proteinsOfInterestNames[proteinIndex], xy = (protein, mRnaExpression[gen][protein] + heighOffset))

		ax.set_title("Generation %s" % gen, fontsize = 12)
		ax.set_xlim([0, numMRnas])
		ax.set_ylim([np.min(mRnaExpression[gen]), np.max(mRnaExpression[gen])])
		ax.spines["top"].set_visible(False)
		ax.tick_params(which = "both", direction = "out", top = "off")
			
	ax.set_xlabel("mRNA transcripts\n(in order of decreasing expected basal expression)", fontsize = 10)

	plt.subplots_adjust(hspace = 1, wspace = 0)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	plt.close("all")

	# Plot for only the proteins of interest
	rows = proteinsOfInterest.shape[0]
	cols = 1
	fig = plt.figure(figsize = (15, 18))

	for proteinIndex, protein in enumerate(proteinsOfInterest):
		ax = plt.subplot(rows, cols, proteinIndex + 1)

		for gen in np.arange(numGens):
			ax.scatter(gen, mRnaExpression[gen][protein])

		ax.set_ylabel(proteinsOfInterestNames[proteinIndex], fontsize = 12)

		ax.spines["top"].set_visible(False)
		ax.tick_params(which = "both", direction = "out", top = "off")	

	ax.set_xlabel("Generation", fontsize = 10)
	plt.subplots_adjust(hspace = 1, wspace = 0)

	exportFigure(plt, plotOutDir, plotOutFileName + "POI", metadata)
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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])