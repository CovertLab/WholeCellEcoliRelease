#!/usr/bin/env python
"""
Plot mRNA half lives (observed vs. actual)

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/30/2015
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

import wholecell.utils.constants
from wholecell.io.tablereader import TableReader

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile))

	isMRna = kb.process.transcription.rnaData["isMRna"]
	isRRna = kb.process.transcription.rnaData["isRRna"]
	isTRna = kb.process.transcription.rnaData["isTRna"]
	rnaIds = kb.process.transcription.rnaData["id"][isMRna]

	expectedDegradationRate = kb.process.transcription.rnaData['degRate'][isMRna].asNumber()


	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	# Note that MoleculeIDs is replaced by objectNames

	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]
	rnaCounts = rnaCountsBulk[1:,:]
	rnaCountsTotal = rnaCounts.sum(axis = 0)

	AllrnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in kb.process.transcription.rnaData["id"]], np.int)
	AllrnaCountsBulk = bulkMolecules.readColumn("counts")[:, AllrnaIndexes]
	AllCounts = AllrnaCountsBulk[1:,:]
	TotalRnaDegraded = (AllCounts * kb.process.transcription.rnaData['degRate'].asNumber()).sum(axis = 1)


	MrnaCounts = AllrnaCountsBulk[1:,isMRna]
	TotalMRnaDegraded = (MrnaCounts * kb.process.transcription.rnaData['degRate'][isMRna].asNumber()).sum(axis = 1)

	RrnaCounts = AllrnaCountsBulk[1:,isRRna]
	TotalRRnaDegraded = (RrnaCounts * kb.process.transcription.rnaData['degRate'][isRRna].asNumber()).sum(axis = 1)

	TrnaCounts = AllrnaCountsBulk[1:,isTRna] 
	TotalTRnaDegraded = (TrnaCounts * kb.process.transcription.rnaData['degRate'][isTRna].asNumber()).sum(axis = 1)


	rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
	time = rnaDegradationListenerFile.readColumn("time")
	countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
	rnaDegradationListenerFile.close()
	rnaDegraded = countRnaDegraded[1:,:]
	rnaDegradedTotal = rnaDegraded.sum(axis = 0)[isMRna]
	rnaDegradationRate = rnaDegradedTotal / 3600. # TODO: this is not true

	rnaSynthesizedListenerFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
	countRnaSynthesized = rnaSynthesizedListenerFile.readColumn('countRnaSynthesized')
	rnaSynthesizedListenerFile.close()
	rnaSynthesized = countRnaSynthesized[1:,:]
	rnaSynthesizedTotal = rnaSynthesized.sum(axis = 0)[isMRna]
	rnaSynthesizedTotalRate = rnaSynthesizedTotal / 3600.
	
	rnaDegradationRate1 = []
	rnaDegradationRate2 = []
	rnaDegradationRate3 = []
	rnaDegradationRate4 = []
	rnaDegradationRate5 = []
	expectedDegradationRateSubset = []
	expectedDegradationRateSubset4 = []
	expectedDegradationRateSubset5 = []

	for i in range(0, len(rnaCountsTotal)): # loop for genes
		if rnaCountsTotal[i] != 0:
			rnaDegradationRate1.append(rnaDegradedTotal[i] / rnaCountsTotal[i]) # Sum_tau(kd*r) / Sum_tau(r)

			rnaDegradationRate2.append(rnaSynthesizedTotal[i] / rnaCountsTotal[i]) # Sum_tau(kb) / Sum_tau(r)

			rnaDegradationRate3.append( (rnaSynthesizedTotal[i] - rnaCounts[0,i]) / rnaCountsTotal[i]) # (Sum_tau(kb) - r) / Sum_tau(r)

			rnaDegradationRate4_t = []
			rnaDegradationRate5_t = []
			for j in range(0, len(rnaDegraded[:,i]) - 1): # loop for timepoints
				if rnaCounts[j,i] != 0:
					if rnaDegraded[j,i] > 0:
						rnaDegradationRate4_t.append(rnaDegraded[j,i] / rnaCounts[j,i]) # Average_t(kd*r / r)
					if (rnaSynthesized[j,i] - (rnaCounts[j+1,i] - rnaCounts[j,i])) > 0:
						rnaDegradationRate5_t.append((rnaSynthesized[j,i] - (rnaCounts[j+1,i] - rnaCounts[j,i])) / rnaCounts[j,i]) # Average_t(kd*r / r)
			
			if len(rnaDegradationRate4_t) > 0:
				rnaDegradationRate4.append(np.mean(rnaDegradationRate4_t))	
				expectedDegradationRateSubset4.append(expectedDegradationRate[i])
			if len(rnaDegradationRate4_t) <= 0:
				rnaDegradationRate4.append(-1)	
				expectedDegradationRateSubset4.append(-1)


			if len(rnaDegradationRate5_t) > 0:
				rnaDegradationRate5.append(np.mean(rnaDegradationRate5_t))
				expectedDegradationRateSubset5.append(expectedDegradationRate[i])	
			if len(rnaDegradationRate5_t) <= 0:
				rnaDegradationRate5.append(-1)	
				expectedDegradationRateSubset5.append(-1)

			expectedDegradationRateSubset.append(expectedDegradationRate[i])

		if rnaCountsTotal[i] == 0:
			rnaDegradationRate1.append(-1)	
			rnaDegradationRate2.append(-1)	
			rnaDegradationRate3.append(-1)	
			expectedDegradationRateSubset.append(-1)
			rnaDegradationRate4.append(-1)	
			expectedDegradationRateSubset4.append(-1)
			rnaDegradationRate5.append(-1)	
			expectedDegradationRateSubset5.append(-1)

	np.savetxt(os.path.join(plotOutDir, 'RNAdecayPredicted.txt'), rnaDegradationRate3)
	np.savetxt(os.path.join(plotOutDir, 'RNAdecayExpected.txt'), expectedDegradationRate)

	# reduction of genes
	expectedDegR = []
	predictedDegR = []
	for i in range(0,len(rnaDegradationRate3)):
		if rnaDegradationRate3[i] != -1 and rnaDegradationRate3[i] != 0 and rnaDegradationRate3[i] != 1:
			if rnaDegradationRate3[i] > 0. and rnaDegradationRate3[i] <= 0.01:
				expectedDegR.append(expectedDegradationRate[i])
				predictedDegR.append(rnaDegradationRate3[i])


	plt.figure(figsize = (8.5, 11))
	maxLine = 1.1 * max(max(expectedDegR), max(predictedDegR))
	
	plt.plot([0, maxLine], [0, maxLine], '--r')	
	plt.plot(expectedDegR, predictedDegR, 'o', markeredgecolor = 'k', markerfacecolor = 'none')
	Correlation_ExpPred = np.corrcoef(expectedDegR, predictedDegR)[0][1]


	plt.xlabel("Expected RNA decay")
	plt.ylabel("Actual RNA decay (at final time step)")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


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
