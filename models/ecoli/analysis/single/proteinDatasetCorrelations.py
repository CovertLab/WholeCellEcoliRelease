#!/usr/bin/env python
"""
Plot WCM protein counts against the data found by Houser et al PLoS CB 2015

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/6/2015
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)


	# Get the names of proteins from the KB
	kb = cPickle.load(open(kbFile, "rb"))

	proteinIds = kb.process.translation.monomerData["id"]

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

	proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]

	bulkMolecules.close()

	trimmedIds = [x[:-3] if x[-3] is '[' else x for x in proteinIds]
	prepedIds = []
	for proteinId in trimmedIds:
		if proteinId == "PROTEIN-CHEA":
			prepedIds.append("null")
		else:
			prepedIds.append(proteinId)

	# Load the Houser Wilke PLoS CB 2015 data from an outside file
	with open("../reconstruction/ecoli/flat/houser_protein_counts.tsv",'r') as f:
		# Dict mapping gene symbol to houser average count
		houserAveCountsDict = {}
		# Read the flat file line by line
		for line in f:
			columns = line.split('\t')
			houserAveCountsDict[columns[0]] = columns[1]

	houserAveCounts = []
	for prepedId in prepedIds:
		if prepedId != "null" and prepedId in kb.moleculeGroups.frameIDGeneSymbol_Dict:
			if kb.moleculeGroups.frameIDGeneSymbol_Dict[prepedId] in houserAveCountsDict:
				houserAveCounts.append(houserAveCountsDict[kb.moleculeGroups.frameIDGeneSymbol_Dict[prepedId]])
				continue
		houserAveCounts.append("null")

	# Load the Taniguichi Xie Science 2010 data from an outside file
	with open("../reconstruction/ecoli/flat/taniguichi_protein_counts.tsv",'r') as f:
		# Dict mapping gene symbol to taniguichi average count
		xieAveCountsDict = {}
		# Read the flat file line by line
		for line in f:
			columns = line.split('\t')
			xieAveCountsDict[columns[0]] = columns[1]

	xieAveCounts = []
	for prepedId in prepedIds:
		if prepedId != "null" and prepedId in kb.moleculeGroups.frameIDGeneSymbol_Dict:
			if kb.moleculeGroups.frameIDGeneSymbol_Dict[prepedId] in xieAveCountsDict:
				xieAveCounts.append(xieAveCountsDict[kb.moleculeGroups.frameIDGeneSymbol_Dict[prepedId]])
				continue
		xieAveCounts.append("null")

	proteinDatasets = [[],[]]

	for idx, houserCount in enumerate(houserAveCounts):
		if houserCount != "null" and xieAveCounts[idx] != "null":
			proteinDatasets[0].append(float(houserCount.strip("\n")))
			proteinDatasets[1].append(float(xieAveCounts[idx].strip("\n")))

	names = ["Houser Wilke", "Taniguichi Xie"]

	plt.figure(figsize= (20,15))

	# Number of whole cell model time points to plot
	timePoints = [100, 1000, 3000]
	for idx, timePoint in enumerate(timePoints):
		wcm_dataset = []
		proteinCounts = proteinCountsBulk[timePoint, :]
		for idx, houserAveCount in enumerate(houserAveCounts):
			if houserAveCount != "null" and xieAveCounts[idx] != "null":
				wcm_dataset.append(proteinCounts[idx])

		proteinDatasets.append(wcm_dataset)		
		names.append("WCM Protein Count (timestep %d)" % timePoint)
	
	from wholecell.analysis.plotting_tools import plotSplom
	plotSplom(proteinDatasets,names)

	
	plt.subplots_adjust(wspace = .7, hspace= .7)
	
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])