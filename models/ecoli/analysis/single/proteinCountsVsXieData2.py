#!/usr/bin/env python
"""
Plot protein monomer counts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
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

# TODO: account for complexation

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

	# Load the Taniguichi Xie Science 2010 data from an outside file
	with open("../reconstruction/ecoli/flat/taniguichi_xie_table_6.tsv",'r') as f:
		# Dict mapping gene symbol to taniguichi average count
		xieAveCountsDict = {}
		# Read the flat file line by line
		for line in f:
			columns = line.split('\t')
			xieAveCountsDict[columns[0]] = columns[4]

	trimmedIds = [x[:-3] if x[-3] is '[' else x for x in proteinIds]
	prepedIds = []
	for proteinId in trimmedIds:
		if proteinId == "PROTEIN-CHEA":
			prepedIds.append("null")
		else:
			prepedIds.append(proteinId)

	xieAveCounts = []
	for prepedId in prepedIds:
		if prepedId != "null" and prepedId in kb.moleculeGroups.frameIDGeneSymbol_Dict:
			if kb.moleculeGroups.frameIDGeneSymbol_Dict[prepedId] in xieAveCountsDict:
				xieAveCounts.append(xieAveCountsDict[kb.moleculeGroups.frameIDGeneSymbol_Dict[prepedId]])
				continue
		xieAveCounts.append("null")

	plt.figure(figsize= (20,15))

	# Number of whole cell model time points to plot
	timePoints = [100, 1000, 3000]
	num_entries = len(timePoints)+1
	proteinCountsMultiTimeStep = []
	names = []

	for idx, timePoint in enumerate(timePoints):

		proteinCounts = proteinCountsBulk[timePoint, :]
		xieAveCountsFinal = []
		proteinCountsFinal = []
		for index, xieAveCount in enumerate(xieAveCounts):
			if xieAveCount != "null":
				xieAveCountsFinal.append(float(xieAveCount))
				proteinCountsFinal.append(proteinCounts[index])
		
		proteinCountsMultiTimeStep.append(proteinCountsFinal)
		names.append("WCM Protein Count (timestep %d)" % timePoint)

	proteinCountsMultiTimeStep.insert(0,xieAveCountsFinal)
	names.insert(0,"Taniguichi Xie Observed Protein Count")

	plottingIndex = 1
	for rowNum in xrange(1,num_entries+1):
		for colNum in xrange(1,num_entries+1):
			if colNum < plottingIndex:
				continue

			corr_coef, pValue = stats.pearsonr(proteinCountsMultiTimeStep[colNum-1], proteinCountsMultiTimeStep[rowNum-1])

			plt.subplot(num_entries,num_entries,num_entries*(rowNum-1)+(colNum))
			plt.scatter(proteinCountsMultiTimeStep[colNum-1], proteinCountsMultiTimeStep[rowNum-1])
			plt.text(0,.9*np.amax(proteinCountsMultiTimeStep[rowNum-1]),"r = %.4f" % (corr_coef), verticalalignment='top', horizontalalignment='left')

			plt.xlabel(names[colNum-1])
			plt.ylabel(names[rowNum-1])

		plottingIndex += 1

	
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
