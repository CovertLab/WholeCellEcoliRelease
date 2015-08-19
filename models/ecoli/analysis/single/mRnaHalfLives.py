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

# import ecocyc_utils

# import tables
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

	# kb = cPickle.load(open(kbFile, "rb"))
	kb = cPickle.load(open(kbFile))

	isMRna = kb.process.transcription.rnaData["isMRna"]
	isRRna = kb.process.transcription.rnaData["isRRna"]
	isTRna = kb.process.transcription.rnaData["isTRna"]
	rnaIds = kb.process.transcription.rnaData["id"][isMRna]
	# import ipdb; ipdb.set_trace()

	expectedDegradationRate = kb.process.transcription.rnaData['degRate'][isMRna].asNumber()

	# counts = RNA(t)
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	# Note that MoleculeIDs is replaced by objectNames
	# import ipdb; ipdb.set_trace()
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]
	rnaCounts = rnaCountsBulk[1:,:]
	rnaCountsTotal = rnaCounts.sum(axis = 0)

	AllrnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in kb.process.transcription.rnaData["id"]], np.int)
	AllrnaCountsBulk = bulkMolecules.readColumn("counts")[:, AllrnaIndexes]
	AllCounts = AllrnaCountsBulk[1:,:]
	TotalRnaDegraded = (AllCounts * kb.process.transcription.rnaData['degRate'].asNumber()).sum(axis = 1)
	# np.savetxt(os.path.join(plotOutDir, 'TotalKdRNA.txt'), TotalRnaDegraded)


	MrnaCounts = AllrnaCountsBulk[1:,isMRna]
	# print kb.process.transcription.rnaData['degRate'][isMRna].asNumber().mean()
	# print kb.process.transcription.rnaData['degRate'][isMRna].asNumber().std()	
	# print len(kb.process.transcription.rnaData['degRate'][isMRna].asNumber())
	# np.savetxt(os.path.join(plotOutDir, 'MRNA.txt'), MrnaCounts.sum(axis = 1))
	TotalMRnaDegraded = (MrnaCounts * kb.process.transcription.rnaData['degRate'][isMRna].asNumber()).sum(axis = 1)
	# np.savetxt(os.path.join(plotOutDir, 'TotalKdMRNA.txt'), TotalMRnaDegraded)

	RrnaCounts = AllrnaCountsBulk[1:,isRRna]
	# print kb.process.transcription.rnaData['degRate'][isRRna].asNumber().mean()
	# print kb.process.transcription.rnaData['degRate'][isRRna].asNumber().std()
	# print len(kb.process.transcription.rnaData['degRate'][isRRna].asNumber())
	# np.savetxt(os.path.join(plotOutDir, 'RRNA.txt'), RrnaCounts.sum(axis = 1))
	TotalRRnaDegraded = (RrnaCounts * kb.process.transcription.rnaData['degRate'][isRRna].asNumber()).sum(axis = 1)
	# np.savetxt(os.path.join(plotOutDir, 'TotalKdRRNA.txt'), TotalRRnaDegraded)

	TrnaCounts = AllrnaCountsBulk[1:,isTRna] 
	# print kb.process.transcription.rnaData['degRate'][isTRna].asNumber().mean()
	# print kb.process.transcription.rnaData['degRate'][isTRna].asNumber().std()
	# print len(kb.process.transcription.rnaData['degRate'][isTRna].asNumber())
	# np.savetxt(os.path.join(plotOutDir, 'TRNA.txt'), TrnaCounts.sum(axis = 1))
	TotalTRnaDegraded = (TrnaCounts * kb.process.transcription.rnaData['degRate'][isTRna].asNumber()).sum(axis = 1)
	# np.savetxt(os.path.join(plotOutDir, 'TotalKdTRNA.txt'), TotalTRnaDegraded)


	# degradation = kd*r(t)
	rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
	time = rnaDegradationListenerFile.readColumn("time")
	countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
	rnaDegradationListenerFile.close()
	rnaDegraded = countRnaDegraded[1:,:]
	rnaDegradedTotal = rnaDegraded.sum(axis = 0)[isMRna]
	rnaDegradationRate = rnaDegradedTotal / 3600. # TODO: this is not true
	# np.savetxt(os.path.join(plotOutDir, 'TotalDegradedMRNA.txt'), rnaDegraded[:, isMRna].sum(axis = 1))
	# np.savetxt(os.path.join(plotOutDir, 'TotalDegradedRRNA.txt'), rnaDegraded[:, isRRna].sum(axis = 1))
	# np.savetxt(os.path.join(plotOutDir, 'TotalDegradedTRNA.txt'), rnaDegraded[:, isTRna].sum(axis = 1))

	# production = kb(t)
	rnaSynthesizedListenerFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
	countRnaSynthesized = rnaSynthesizedListenerFile.readColumn('countRnaSynthesized')
	rnaSynthesizedListenerFile.close()
	rnaSynthesized = countRnaSynthesized[1:,:]
	rnaSynthesizedTotal = rnaSynthesized.sum(axis = 0)[isMRna]
	rnaSynthesizedTotalRate = rnaSynthesizedTotal / 3600.
	# np.savetxt(os.path.join(plotOutDir, 'TotalSynthesizedMRNA.txt'), rnaSynthesized[:, isMRna].sum(axis = 1))
	# np.savetxt(os.path.join(plotOutDir, 'TotalSynthesizedRRNA.txt'), rnaSynthesized[:, isRRna].sum(axis = 1))
	# np.savetxt(os.path.join(plotOutDir, 'TotalSynthesizedTRNA.txt'), rnaSynthesized[:, isTRna].sum(axis = 1))
	
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

	# np.savetxt(os.path.join(plotOutDir, 'RNAdecayPredicted1.txt'), rnaDegradationRate1)
	# np.savetxt(os.path.join(plotOutDir, 'RNAdecayPredicted2.txt'), rnaDegradationRate2)
	np.savetxt(os.path.join(plotOutDir, 'RNAdecayPredicted3.txt'), rnaDegradationRate3)
	# np.savetxt(os.path.join(plotOutDir, 'RNAdecayPredicted4.txt'), rnaDegradationRate4)
	# np.savetxt(os.path.join(plotOutDir, 'RNAdecayPredicted5.txt'), rnaDegradationRate5)
	np.savetxt(os.path.join(plotOutDir, 'RNAdecayExpected.txt'), expectedDegradationRate)

	# print np.corrcoef(expectedDegradationRateSubset, rnaDegradationRate1)[0][1]
	# print np.corrcoef(expectedDegradationRateSubset, rnaDegradationRate2)[0][1]
	# print np.corrcoef(expectedDegradationRateSubset, rnaDegradationRate3)[0][1]
	# print np.corrcoef(expectedDegradationRateSubset4, rnaDegradationRate4)[0][1]
	# print np.corrcoef(expectedDegradationRateSubset5, rnaDegradationRate5)[0][1]

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
	#plt.errorbar(expectedDegradationRate, rnaDegradationRate, rnaDegradationRateStd)
	Correlation_ExpPred = np.corrcoef(expectedDegR, predictedDegR)[0][1]

	# Saving data files for predicted and observed RNA decays and B-numbers for all mRNAs
	# bNumbersToFrameIds = ecocyc_utils.get_bNumbersToFrameIds("frame-id-names.txt")
	# frameIdsToBNumbers = ecocyc_utils.get_frameIdsToBNumbers("frame-id-names.txt")
	# mRNAgeneIds = kb.process.transcription.rnaData["geneId"][isMRna]
	# mRNAgeneBnumber = []
	# for i in range(0,len(mRNAgeneIds)):
	# 	mRNAgeneBnumber = np.append(mRNAgeneBnumber, frameIdsToBNumbers[mRNAgeneIds[i]])
	# text_file = open(os.path.join(plotOutDir, 'mRNABnumber.txt'), "w")
	# text_file.writelines(["%s\n" % item  for item in mRNAgeneBnumber])
	# text_file.close()

	print "Correlation expected and predicted half-lives = %.3f" % Correlation_ExpPred

	LifetimePred_avg = 1 / np.average(rnaDegradationRate) / 60
	LifetimePred_std = 1 / np.std(rnaDegradationRate) / 60
	# print "Avg predicted lifetime = %.3f min" % LifetimePred_avg
	# print "Std predicted lifetime = %.3f min" % LifetimePred_std

	LifetimeExp_avg = 1 / np.average(expectedDegradationRate) / 60
	LifetimeExp_std = 1 / np.std(expectedDegradationRate) / 60
	# print "Avg expected lifetime = %.3f min" % LifetimeExp_avg
	# print "Std expected lifetime = %.3f min" % LifetimeExp_std
	
	ResSynthRate_avg = np.average(rnaDegradationRate - expectedDegradationRate)
	# print ResSynthRate_avg
	ResSynthRate_std = np.std(rnaDegradationRate - expectedDegradationRate)
	# print ResSynthRate_std
	LifetimeRes_avg = 1 / ResSynthRate_avg / 60
	LifetimeRes_std = 1 / ResSynthRate_std / 60
	# print "Avg residual lifetime = %.3f min" % LifetimeRes_avg 
	# print "Std residual lifetime = %.3f min" % LifetimeRes_std

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
